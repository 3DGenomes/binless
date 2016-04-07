////
// GAM : N ~ s(log(distance)) + is.close + log(rejoined.1) + log(rejoined.2)
// using cubic splines and 2nd order penalty with log-link and negative binomial
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

  // design matrix for count ~ s(log(dist)) + isclose + log(rejoined.1) + log(rejoined.2)
  // X is N x (K+3)
  matrix design(vector dist, int[] isclose, vector rejoined1, vector rejoined2, int K, int q) {
    ////intercept, centered bspline basis with K-1 params, one binary factor and two regressors
    matrix[rows(dist),K] X;
    X <- bspline(log(dist), K-q, q);
    return append_col(rep_vector(1,rows(dist)),
                      append_col(center(X,X),
                                 append_col(to_vector(isclose),
                                            append_col(log(rejoined1),log(rejoined2)))));
  }
  
  // difference operator for that design matrix
  matrix difference_op(vector dist, int K, int q, int d) {
    matrix[K-d,K-1] diff;
    diff <- center(bspline(log(dist), K-q, q), difference_matrix_sqrt(K,d));
    return append_col(rep_vector(0,K-d), 
                      append_col(diff, rep_matrix(0,K-d,3))); //intercept, factor and regressors have no constraint
  }
  
}
////////////////////////////////////////////////////////////////
data {
  int<lower=1> N; //number of data points
  int<lower=1> K; //number of parameters
  int count[N]; //dependent variable
  vector[N] dist; //genomic distance for this count
  vector[N] rejoined1; //independent variables
  vector[N] rejoined2;
  int<lower=0,upper=1> isclose[N];
}
transformed data {
  matrix[N,K+3] X; //design matrix
  matrix[K-difforder(),K+3] P;
  X <- design(dist, isclose, rejoined1, rejoined2, K, splinedegree());
  P <- difference_op(dist, K, splinedegree(), difforder());
}
parameters {
  vector[K+3] beta;
  real<lower=0> alpha;
  real<lower=0> lambda;
}
model {
  //negative binomial GAM with log link
  count ~ neg_binomial_2_log(X * beta, alpha);
  //P-spline prior on the differences (K-1 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  P*beta ~ normal(0, 1./(alpha*lambda));
}
generated quantities {
  matrix[N,K+3] designmat; //design matrix
  matrix[N,K+3] weighted; //weighted basis functions
  vector[N] pred; //mean predictive
  vector[N] decay; //spline interpolant
  vector[K+3] edfvec;
  real edf; //effective degrees of freedom
  designmat <- X;
  weighted <- X .* rep_matrix(beta', rows(X));
  pred <- exp(X * beta);
  decay <- exp(X[,2:K]*beta[2:K]);
  {
    matrix[K+3,K+3] XtX;
    matrix[K+3,K+3] jitter;
    XtX <- crossprod(X);
    jitter <- diag_matrix(rep_vector(mean(diagonal(XtX))*0.01,K+3));
    edfvec <- diagonal(inverse_spd(jitter+XtX+lambda*crossprod(P)) * XtX);
    edf <- sum(edfvec);
  }
}
