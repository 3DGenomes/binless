////
// GAM : y ~ s(x1)+s(x2) using cubic splines
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
    
  //apply centering constraint on a design matrix
  matrix center(matrix X) {
    int N;
    int K;
    N <- rows(X);
    K <- cols(X);
    {
      //form QR decomposition of Xt1
      matrix[K,1] sums;
      vector[K] u;
      sums <- to_matrix(X'*rep_vector(1,N));
      u <- householder(sums)[,1];
      return X[,2:K] - (2*X*u)*u[2:K]';
    }
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

  // design matrix for 1D b-spline fits
  // K1 params for s(x1), K2 params for s(x2), q is spline degree
  // design matrix is N x (K1+K2-1) e.g. intercept and centering constraints
  matrix design(vector x1, vector x2, int K1, int K2, int q) {
    return append_col(rep_vector(1,rows(x1)), 
                      append_col(center(bspline(x1, K1-q, q)), center(bspline(x2, K2-q, q))));
  }
  
  // difference operator for that design matrix
  // K1 params for s(x1), K2 params for s(x2), d is difference degree
  // difference operator is AxB where
  //   A is sum_i Ki-1-d
  //   B is 1 + sum_i Ki-1
  matrix difference_op(int K1, int K2, int d) {
    matrix[K1+K2-(d+1)*2,K1+K2-1] ret;
    ret <- rep_matrix(0,rows(ret),cols(ret));
    ret[1:(K1-1-d),2:K1] <- difference_matrix_sqrt(K1-1,d);
    ret[(K1-d):,(K1+1):] <- difference_matrix_sqrt(K2-1,d);
    return ret; 
  }
  
  //design matrix for penalty
  // K1 params for s(x1), K2 params for s(x2), d is difference degree
  // design matrix is Ax2 where
  //   A is sum_i Ki-1-d
  matrix design_penalty(int K1, int K2, int d) {
    matrix[K1+K2-(d+1)*2,2] ret;
    ret <- rep_matrix(0,K1+K2-(d+1)*2,2);
    ret[:(K1-1-d),1] <- rep_vector(1,K1-1-d);
    ret[(K1-d):,2] <- rep_vector(1,K2-1-d);
    return(ret);
  }
}
////////////////////////////////////////////////////////////////
data {
  int<lower=1> N; //number of data points
  int<lower=1> K1; //number of parameters for s1
  int<lower=1> K2;
  int y[N]; //dependent variable
  vector[N] x1; //independent variables
  vector[N] x2;
}
transformed data {
  matrix[N,K1+K2-1] X; //design matrix
  matrix[K1+K2-(difforder()+1)*2, K1+K2-1] P; //difference matrix
  matrix[K1+K2-(difforder()+1)*2, 2] Z; //design matrix for penalty part
  X <- design(x1, x2, K1, K2, splinedegree());
  P <- difference_op(K1,K2,difforder());
  Z <- design_penalty(K1,K2,difforder());
}
parameters {
  vector[K1+K2-1] beta;
  real<lower=0> alpha;
  vector<lower=0>[2] lambda;
}
model {
  //exponential GAM
  y ~ neg_binomial_2(X * beta, alpha);
  //P-spline prior on the differences
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  P*beta ~ normal(0, 1. ./ (Z*lambda*alpha));
}
generated quantities {
  matrix[N,K1+K2-1] weighted; //weighted basis functions
  vector[N] pred; //predicted mean
  real intercept;
  vector[N] splines[2]; //fits of each component
  vector[K1+K2-1] edfvec; //effective degrees of freedom per parameter
  real edf; //effective degrees of freedom
  //
  weighted <- X .* rep_matrix(beta', rows(X));
  pred <- X * beta;
  intercept <- beta[1];
  splines[1] <- X[,2:K1]*beta[2:K1];
  splines[2] <- X[,(K1+1):]*beta[(K1+1):];
  {
    //dof matrix is ((Xt*X + Pt * Z*lambda *P)^-1 * Xt*X)
    //assumes likelihood is separable normal
    matrix[K1+K2-1,K1+K2-1] XtX;
    matrix[K1+K2-1,K1+K2-1] jitter;
    XtX <- crossprod(X);
    jitter <- diag_matrix(rep_vector(mean(diagonal(XtX))*0.01,K1+K2-1));
    edfvec <- diagonal(inverse_spd(jitter+XtX+quad_form_sym(diag_matrix(Z*lambda), P)) * XtX);
    edf <- sum(edfvec);
  }
}
