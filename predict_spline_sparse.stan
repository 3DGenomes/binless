////
// predict a sparse cubic spline GAM model: y ~ s(x) (with intercept)
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
  int K; //number of functions in spline base
  int<lower=1> N; //number of samples along x
  vector[2] xrange; //first and last values of x
  //fitted things
  real intercept;
  vector[K-1] beta;
}
transformed data {
  vector[N] gen; //x positions, has to be within bounds of xrange
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(N)] Xrow_w;
  int Xrow_v[nnz(N)];
  int Xrow_u[N+1];
  vector[K] prow;
  
  ////bias spline, sparse
  gen <- xrange[1] + range(0,N-1)*(xrange[2]-xrange[1])/(N-1);
  //BEGIN sparse calculation
  //cannot write function that modifies its arguments so we put it here
  //input: vector[N] genome, int K, int splinedegree()
  //output: vector[nnz(N)] Xrow_w, int Xrow_v[nnz(N)], int Xrow_u[N+1]
  {
    real dx; //interval length
    row_vector[K] t; //Knot locations (except last)
    int x_n[K-splinedegree()+1]; //cumulative histogram of genome values
    int idx_u; //counters for filling of sparse matrix
    int idx_w;
    //
    dx <- 1.01*(max(xrange)-min(xrange))/(K-splinedegree()); //maKe it slightly larger
    t <- min(xrange) - dx*0.01 + dx * range(-splinedegree(),K-splinedegree()-1)';
    //get indices of genome which cross to the next segment and build x_n.
    x_n[1] <- 1;
    x_n[2:(K-splinedegree())] <- cumulative_hist(gen, t[(splinedegree()+2):]);
    x_n[K-splinedegree()+1] <- rows(gen)+1;
    //
    //build spline, interval per interval
    idx_u <- 1;
    idx_w <- 1;
    for (i in 1:(K-splinedegree())) {
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
  
  //projector for biases (GAM)
  prow <- column_sums(K, Xrow_w, Xrow_v)';
  prow <- prow / sqrt(dot_self(prow));
}
parameters {}
model {}
generated quantities {
  vector[N] x;
  vector[N] log_mean; //log of mean parameter
  vector[K] beta_centered; //weights to the unconstrained spline base
  matrix[N,4] basis; //basis functions, wrapped
  matrix[N,4] weighted; //same but weighted
  
  x <- gen;
  //mean
  {
    vector[K] beta_aug;
    beta_aug[1] <- sum(beta);
    beta_aug[2:] <- beta;
    beta_centered <- beta_aug - (beta_aug' * prow) * prow;
    log_mean <- csr_matrix_times_vector(N, K, Xrow_w, Xrow_v, Xrow_u, beta_centered);
  }
  //basis functions
  {
    vector[K] bt;
    vector[4] pattern;
    int patternsz;
    int npatterns;
    patternsz <- 4; //splinedegree()+1;
    npatterns <- K/patternsz; // K = npatterns*patternsz + remainder
    //first pattern
    pattern[1] <- 1;
    pattern[2:patternsz] <- rep_vector(0,patternsz-1);
    //inner patterns
    for (i in 1:npatterns) {
      bt[((i-1)*patternsz+1):(i*patternsz)] <- pattern;
    }
    //last pattern
    if (K-npatterns*patternsz>0) bt[npatterns*patternsz:] <- pattern[1:(K-npatterns*patternsz)];
    //generate wrapped bases
    for (i in 1:patternsz) {
      vector[K] tmp;
      basis[,i] <- csr_matrix_times_vector(N, K, Xrow_w, Xrow_v, Xrow_u, bt);
      weighted[,i] <- csr_matrix_times_vector(N, K, Xrow_w, Xrow_v, Xrow_u, bt .* beta_centered);
      tmp[1] <- bt[K];
      tmp[2:] <- bt[:(K-1)];
      bt <- tmp;
    }
  }
}