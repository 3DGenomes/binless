////
// GAM : y ~ s(X) using cubic splines with log-link
////
functions {
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
  }
  int splinedegree() {return 3;} //set 3 for cubic spline
  int difforder() {return 2;} //set 2 for 2nd order difference penalty

  // Evaluate b-spline basis functions of degree q at the values
  // in x within support [xmin,xmax] and with ndx intervals.
  // Use q=3 for cubic spline.
  // adapted from Eilers and Marx, Stat Sci, 1996
  matrix bspline(vector x, int ndx, int q) {
    real dx;
    row_vector[ndx+q] t;
    int r[ndx+q];
    matrix[rows(x), ndx+q] T;
    matrix[rows(x), ndx+q] X;
    matrix[rows(x), ndx+q] P;
    matrix[rows(x), ndx+q] B;
    dx <- 1.01*(max(x)-min(x))/ndx; //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,ndx-1)';
    for (i in 2:ndx+q) r[i-1] <- i;
    r[ndx+q] <- 1;
    T <- rep_matrix(t, rows(x));
    X <- rep_matrix(x, ndx+q);
    P <- (X - T) / dx;
    for (i in 1:rows(x))
      for (j in 1:cols(t))
        B[i,j] <- (T[i,j] <= X[i,j]) && (X[i,j] < T[i,j]+dx);
    for (k in 1:q)
      B <- ( P .* B + (k+1-P) .* B[,r]) / k;
    return B;
  }
  
  //collection of householder vectors to decompose A
  //if A is n x m then the function returns a n x 2m matrix [U : R]
  // U will be n x m and the vectors are stored in the lower triangle of U.
  // Q is prod_i 1 - 2u_iu_iT 
  // R is upper triangular.
  matrix householder(matrix A) {
    int n;
    int m;
    n <- rows(A);
    m <- cols(A);
    {
      matrix[n,m] U;
      matrix[n,m] R;
      vector[n] e;
      U <- rep_matrix(0,n,m);
      e[2:n] <- rep_vector(0,n-1);
      e[1] <- 1;
      R <- A;
      for (k in 1:m) {
        vector[n-k+1] x;
        vector[n-k+1] u;
        x <- R[k:n,k];
        u <- sqrt(x'*x)*e[1:(n-k+1)] + x;
        if (x[1]<0) u <- -u;
        u <- u/sqrt(u'*u);
        {
          matrix[n-k+1,m-k+1] tmp; //stan 2.9.0 issues compile error for deep_copy
          tmp <- R[k:n,k:m] - 2*u*transpose(u)*R[k:n,k:m];
          R[k:n,k:m] <- tmp;
        }
        U[k:n,k] <- u;
      }
      return append_col(U,R);
    }
  }
    
  //compute householder vector of centering constraint X^T * 1
  vector centering_constraint(matrix X) {
    int N;
    int K;
    N <- rows(X);
    K <- cols(X);
    {
      //form QR decomposition of Xt1
      matrix[K,1] sums;
      vector[K] u;
      sums <- to_matrix(rep_row_vector(1,N)*X)';
      return householder(sums)[,1];
    }
  }
  
  //apply centering constraint on X to a matrix D: D -> DZ
  // where Z is such that X^T 1 = QR with Q=[u:Z]
  // if D is LxK, and X is NxK, returns a Lx(K-1) matrix
  matrix center(matrix X, matrix D) {
    vector[cols(X)] u;
    u <- centering_constraint(X);
    return D[,2:] - (2*D*u)*u[2:]';
  }
  
  //difference penalty for P-splines
  //K is the size of the parameter vector
  //d is the order of the difference
  //Returns P, which is K-d x K
  //D=t(P)P is the difference matrix, P*beta is the difference vector
  matrix difference_matrix_sqrt(int K, int d) {
    matrix[K,K] P;
    P <- diag_matrix(rep_vector(1, K));
    for (i in 1:d) {
      matrix[K-i,K] tmp; //stan 2.9.0 issues compile error for deep_copy
      tmp <- P[2:(K-i+1),]-P[1:(K-i),];
      P[1:(K-i),] <- tmp;
    }
    return P[1:(K-d),];
  }

  // design matrix for 1D b-spline fit of degree q with K functions
  matrix design(vector x, int K, int q) {
    matrix[rows(x),K] X;
    X <- bspline(x, K-q, q);
    return append_col(rep_vector(1,rows(x)), X);
  }
  
  // difference operator for that design matrix
  matrix difference_op(vector x, int K, int q, int d) {
    return difference_matrix_sqrt(K,d);
  }
}
////////////////////////////////////////////////////////////////
data {
  int<lower=1> N; //number of data points
  int<lower=1> K; //number of parameters
  vector[N] y; //dependent variable
  vector[N] x; //independent variable
}
transformed data {
  matrix[N,K+1] X; //design matrix
  matrix[K-difforder(),K] P;
  vector[K] spline_hh;
  X <- design(x, K, splinedegree());
  P <- difference_op(x, K, splinedegree(), difforder());
  {
    vector[K] spline_v;
    spline_v <- centering_constraint(X[,2:(K+1)]);
    spline_hh <- -2*spline_v[1]*spline_v;
    spline_hh[1] <- 1+spline_hh[1];
  }
}
parameters {
  real intercept;
  ordered[K] beta_spline;
  real<lower=0> sigma2;
  real<lower=0> lambda;
}
transformed parameters {
  vector[K+1] beta; //intercept and centered spline
  beta[1] <- intercept;
  beta[2:] <- beta_spline - (spline_hh'*beta_spline)*spline_hh; //center the coefs
}
model {
  //exponential GAM
  y ~ normal(X * beta, sigma2);
  //P-spline prior on the differences (K-1 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  P*beta_spline ~ normal(0, sigma2/lambda);
}
generated quantities {
  matrix[N,K+1] designmat; //design matrix
  matrix[N,K+1] weighted; //weighted basis functions
  vector[N] pred; //spline interpolant
  real edf; //effective degrees of freedom
  vector[K+1] edfvec;
  designmat <- X;
  weighted <- X .* rep_matrix(beta', rows(x));
  pred <- (X * beta);
  {
    matrix[K+1,K+1] XtX;
    matrix[K+1,K+1] S;
    XtX <- crossprod(X);
    S[2:,2:] <- crossprod(P);
    S[1,1] <- 1;
    S[2:,1] <- rep_vector(0,K);
    S[1,2:] <- rep_row_vector(0,K);
    edfvec <- diagonal(inverse_spd(XtX+lambda*S) * XtX);
    edf <- sum(edfvec);
  }
}
