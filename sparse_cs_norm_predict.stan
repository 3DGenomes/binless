////
// Cut-site normalization model: predict expected biases/counts at new locations
////
functions {
  ///////////
  //BEGIN internal use
  //
  //generate [imin, imin+1, ... imax]
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
  }
  //Compute cumulative histogram of x at left cut points q
  //i.e. index of first x value that falls in bins defined by cutpoints defined in q
  //assumes x is ordered.
  //x[hist[i]] is the first value in x that is >= q[i]
  int[] cumulative_hist(vector x, row_vector q) {
    int indices[cols(q)];
    int ix;
    int N;
    ix <- 1;
    N <- rows(x);
    for (iq in 1:cols(q)) {
      while (ix < N && x[ix] < q[iq]) ix <- ix + 1;
      indices[iq] <- ix;
    }
    return indices;
  }
  matrix bspl_gen(vector x, real dx, row_vector t, int q) {
    int N;
    int K;
    K <- cols(t);
    N <- rows(x);
    {
      int r[K];
      matrix[N, K] T;
      matrix[N, K] X;
      matrix[N, K] P;
      matrix[N, K] B;
      for (i in 2:K) r[i-1] <- i;
      r[K] <- 1;
      T <- rep_matrix(t, N);
      X <- rep_matrix(x, K);
      P <- (X - T) / dx;
      for (i in 1:N)
        for (j in 1:K)
          B[i,j] <- (T[i,j] <= X[i,j]) && (X[i,j] < T[i,j]+dx);
      for (k in 1:q)
        B <- ( P .* B + (k+1-P) .* B[,r]) / k;
      return B;
    }
  }
  //END internal use
  ///////////
  
  int splinedegree() {return 3;} //set 3 for cubic spline
  
  int nnz(int N) {return N*(splinedegree()+1);} //nonzero count for design matrix
  
  row_vector column_sums(int K, vector Xw, int[] Xv) {
    row_vector[K] sums;
    sums <- rep_row_vector(0,K);
    for (i in 1:rows(Xw)) {
      sums[Xv[i]] <- sums[Xv[i]] + Xw[i];
    }
    return sums;
  }
  
  matrix bspline(vector x, int K, int q) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(max(x)-min(x))/(K-q); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }


}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Krow; //number of functions in spline base for row biases
  int Kdiag; //number of functions in spline base for diagonal decay
  //biases
  vector[2] cutsites; //first and last cut site locations, used for range
  int<lower=1> S; //number of samples along genome
  //counts
  int<lower=1> N; //number of samples along diagonal
  //fitted things
  real intercept;
  real eRJ; //exposure for rejoined ends
  real eDE; //exposure for dangling ends
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  positive_ordered[Kdiag-1] beta_diag;
}
transformed data {
  vector[S] gen; //where to look at on the genome, has to be within bounds of cutsites
  vector[N] dist; //where to look at in the decay, has to be <= cutsites[2]-cutsites[1]
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[Krow] prow;
  //diagonal SCAM spline, dense
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag-1] pdiag;
  
  ////bias spline, sparse (nu and delta have the same design)
  gen <- cutsites[1] + range(0,S-1)*(cutsites[2]-cutsites[1])/(S-1);
  //BEGIN sparse calculation
  //cannot write function that modifies its arguments so we put it here
  //input: vector[N] genome, int Krow, int splinedegree()
  //output: vector[nnz(S)] Xrow_w, int Xrow_v[nnz(S)], int Xrow_u[S+1]
  {
    real dx; //interval length
    row_vector[Krow] t; //Krownot locations (except last)
    int x_n[Krow-splinedegree()+1]; //cumulative histogram of genome values
    int idx_u; //counters for filling of sparse matrix
    int idx_w;
    //
    dx <- 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree()); //maKrowe it slightly larger
    t <- min(cutsites) - dx*0.01 + dx * range(-splinedegree(),Krow-splinedegree()-1)';
    //get indices of genome which cross to the next segment and build x_n.
    x_n[1] <- 1;
    x_n[2:(Krow-splinedegree())] <- cumulative_hist(gen, t[(splinedegree()+2):]);
    x_n[Krow-splinedegree()+1] <- rows(gen)+1;
    //
    //build spline, interval per interval
    idx_u <- 1;
    idx_w <- 1;
    for (i in 1:(Krow-splinedegree())) {
      //at any genome, there are splinedegree()+1 nonzero b-splines. Compute them.
      int xbegin;
      int xend;
      xbegin <- x_n[i];
      xend <- x_n[i+1]-1;
      {
        matrix[xend-xbegin+1,splinedegree()+1] tmp;
        tmp <- bspl_gen(gen[xbegin:xend], dx, t[i:(i+splinedegree())], splinedegree());
        for (ti in 1:(xend-xbegin+1)) {
          Xrow_u[idx_u] <- idx_w;
          idx_u <- idx_u + 1;
          for (tj in 1:(splinedegree()+1)) {
            Xrow_w[idx_w] <- tmp[ti,tj];
            Xrow_v[idx_w] <- i+tj-1;
            idx_w <- idx_w + 1;
          }
        }
      }
    }
    Xrow_u[idx_u] <- idx_w; 
  }
  //END sparse calculation

  //diagonal SCAM spline, dense
  dist <- range(1,N)*(cutsites[2]-cutsites[1])/N; //avoid 0
  {
    Xdiag <- bspline(log(dist), Kdiag, splinedegree());
  }
  
  
  //projector for biases (GAM)
  prow <- column_sums(Krow, Xrow_w, Xrow_v)';
  prow <- prow / sqrt(dot_self(prow));
  //projector for diagonal (SCAM)
  {
    row_vector[Kdiag] tmp;
    tmp <- rep_row_vector(1,N) * Xdiag;
    pdiag <- -tmp[2:] / (tmp * rep_vector(1,Kdiag));
  }
}
parameters {}
model {}
generated quantities {
  //viewpoints
  vector[S] genome;
  vector[N] distnce;
  //nu
  vector[S] log_nu; // log(nu)
  //delta
  vector[S] log_delta; // log(delta)
  //diag
  vector[N] log_decay;
  //means
  vector[S] mean_DL;
  vector[S] mean_DR;
  vector[S] mean_RJ;
  
  genome <- gen;
  distnce <- dist;
  //nu
  {
    vector[Krow] beta_nu_centered;
    vector[Krow] beta_nu_aug;
    beta_nu_aug[1] <- sum(beta_nu);
    beta_nu_aug[2:] <- beta_nu;
    beta_nu_centered <- beta_nu_aug - (beta_nu_aug' * prow) * prow;
    log_nu <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
  }
  //delta
  {
    vector[Krow] beta_delta_centered;
    vector[Krow] beta_delta_aug;
    beta_delta_aug[1] <- sum(beta_delta);
    beta_delta_aug[2:] <- beta_delta;
    beta_delta_centered <- beta_delta_aug - (beta_delta_aug' * prow) * prow;
    log_delta <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
  }
  //decay
  {
    vector[Kdiag] beta_diag_centered;
    real epsilon;
    real val;
    epsilon <- -1; //+1 for increasing, -1 for decreasing spline
    val <- epsilon*pdiag*beta_diag;
    beta_diag_centered[1] <- val;
    beta_diag_centered[2:] <- val+epsilon*beta_diag;
    log_decay <- Xdiag * beta_diag_centered;
  }
   //means
  {
    vector[S] base_bias;
    vector[N] base_count;
    //
    base_bias <- intercept + log_nu;
    mean_RJ <- exp(base_bias + eRJ);
    mean_DL <- exp(base_bias + eDE + log_delta);
    mean_DR <- exp(base_bias + eDE - log_delta);
  }
}