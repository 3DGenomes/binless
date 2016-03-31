////
// GAM : y ~ s(X) using cubic splines with log-link
// The design matrix is the original spline base, sparse storage
// Intercept is modelled explicitly, and the spline is centered at each
// iteration. Secondary order differences are used.
// EDF calculation is not performed because no sparse Cholesky exists in stan.
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
}
////////////////////////////////////////////////////////////////
data {
  int<lower=1> N; //number of data points
  int<lower=1> K; //number of parameters
  int y[N]; //dependent variable
  vector[N] x; //independent variable
}
transformed data {
  vector[nnz(N)] Xw; //sparse design matrix for spline
  int Xv[nnz(N)];
  int Xu[N+1];
  vector[K] p; //projector on centering constraint
  //BEGIN sparse calculation
  //cannot write function that modifies its arguments so we put it here
  //input: vector[N] x, int K, int q
  //output: vector[nnz(N)] Xw, int Xv[nnz(N)], int Xu[N+1]
  {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    int x_n[K-splinedegree()+1]; //cumulative histogram of x values
    int idx_u; //counters for filling of sparse matrix
    int idx_w;
    //
    dx <- 1.01*(max(x)-min(x))/(K-splinedegree()); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-splinedegree(),K-splinedegree()-1)';
    //get indices of x which cross to the next segment and build x_n.
    x_n[1] <- 1;
    x_n[2:(K-splinedegree())] <- cumulative_hist(x, t[(splinedegree()+2):]);
    x_n[K-splinedegree()+1] <- rows(x)+1;
    //
    //build spline, interval per interval
    idx_u <- 1;
    idx_w <- 1;
    for (i in 1:(K-splinedegree())) {
      //at any x, there are splinedegree()+1 nonzero b-splines. Compute them.
      int xbegin;
      int xend;
      xbegin <- x_n[i];
      xend <- x_n[i+1]-1;
      {
        matrix[xend-xbegin+1,splinedegree()+1] tmp;
        tmp <- bspl_gen(x[xbegin:xend], dx, t[i:(i+splinedegree())], splinedegree());
        for (ti in 1:(xend-xbegin+1)) {
          Xu[idx_u] <- idx_w;
          idx_u <- idx_u + 1;
          for (tj in 1:(splinedegree()+1)) {
            Xw[idx_w] <- tmp[ti,tj];
            Xv[idx_w] <- i+tj-1;
            idx_w <- idx_w + 1;
          }
        }
      }
    }
    Xu[idx_u] <- idx_w; 
  }
  //END sparse calculation
  //build projector on centering constraint
  p <- column_sums(K, Xw, Xv)';
  p <- p / sqrt(dot_self(p));
}
parameters {
  real intercept;
  vector[K-1] beta;
  real<lower=0> alpha;
  real<lower=0> lambda;
}
transformed parameters {
  vector[K] beta_centered;
  vector[K] beta_aug;
  beta_aug[1] <- sum(beta);
  beta_aug[2:] <- beta;
  beta_centered <- beta_aug - (beta_aug' * p) * p;
}
model {
  //exponential GAM
  y ~ neg_binomial_2(exp(intercept + csr_matrix_times_vector(N, K, Xw, Xv, Xu, beta_centered)), alpha);
  //P-spline prior on the 2nd order differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_aug[:(K-2)]-2*beta_aug[2:(K-1)]+beta_aug[3:] ~ normal(0, 1./(alpha*lambda));
}
generated quantities {
  real offset;
  vector[N] pred; //spline interpolant
  offset <- intercept;
  pred <- exp(intercept + csr_matrix_times_vector(N, K, Xw, Xv, Xu, beta_centered));
}
