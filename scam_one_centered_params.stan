////
// SCAM : y ~ s(X) using monotone decreasing cubic splines with log-link
// The design matrix is the original spline base,
// Intercept is modelled explicitly, and the spline is centered at each
// iteration
////
functions {
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
  }
  int splinedegree() {return 3;} //set 3 for cubic spline
  int difforder() {return 2;} //set 2 for 2nd order difference penalty

  // Evaluate b-spline basis functions of degree q at the values
  // in x within support [xmin,xmax] and with K basis functions.
  // Use q=3 for cubic spline.
  // adapted from Eilers and Marx, Stat Sci, 1996
  matrix bspline(vector x, int K, int q) {
    real dx;
    row_vector[K] t;
    int r[K];
    matrix[rows(x), K] T;
    matrix[rows(x), K] X;
    matrix[rows(x), K] P;
    matrix[rows(x), K] B;
    dx <- 1.01*(max(x)-min(x))/(K-q); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,K-q-1)';
    for (i in 2:K) r[i-1] <- i;
    r[K] <- 1;
    T <- rep_matrix(t, rows(x));
    X <- rep_matrix(x, K);
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
  
  //Sigma matrix for SCAM smoothing spline with K params.
  //Sigma is K x K
  //t=1: monotone increasing
  //t=2: monotone decreasing
  matrix scam_sigma(int K, int t);
  matrix scam_sigma(int K, int t) {
    matrix[K,K] Sigma;
    if (t==1) { //monotone increasing
      Sigma[,1] <- rep_vector(1,K);
      for (j in 2:K) {
        Sigma[:(j-1),j] <- rep_vector(0,j-1);
        Sigma[j:,j] <- rep_vector(1,K-j+1);
      }
    } else if (t==2) { //monotone decreasing
      Sigma[2:K,2:K] <- -scam_sigma(K-1, 1);
      Sigma[,1] <- rep_vector(1,K);
      Sigma[1,2:K] <- rep_row_vector(0,K-1);
    }
    return Sigma;
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
  //difference penalty for SCAM
  //K is the size of the parameter vector
  //Returns P, which is K-2 x K
  //D=t(P)P is the difference matrix, P*beta is the difference vector
  //t=1: monotone increasing
  //t=2: monotone decreasing
  matrix scam_difference_op(int K, int t) {
    matrix[K-2,K] P;
    P <- rep_matrix(0,K-2,K);
    if (t == 1 || t == 2) {
      for (i in 1:(K-2)){
        P[i,i+1] <- 1;
        P[i,i+2] <- -1;
      }
    }
    return P;
  }


  // design matrix for 1D b-spline fit of degree q with K functions
  matrix design_spline(vector x, int K, int q) {
    return bspline(x, K, q);
  }
  
  // difference operator for that design matrix
  matrix difference_op(int K, int d) {
    return scam_difference_op(K,d); //difference_matrix_sqrt(K,d);
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
  matrix[N,K] Xs; //design matrix for spline
  matrix[K-difforder(),K] Diff; //difference operator for penalty
  row_vector[K-1] p; //projector on centering constraint
  matrix[K,K] Sigma;
  Xs <- design_spline(x, K, splinedegree());
  Diff <- difference_op(K, 2); //difforder());
  Sigma <- scam_sigma(K, 2); //monotone decreasing
  {
    row_vector[K] tmp;
    tmp <- rep_row_vector(1,N) * Xs;
    p <- (-tmp * Sigma[,2:]) / (tmp * Sigma[,1]);
  }
}
parameters {
  real intercept;
  vector[K-1] beta;
  real<lower=0> sigma2;
  real<lower=0> lambda;
}
transformed parameters {
  vector[K] beta_aug;
  vector[K] beta_centered;
  beta_aug[2:] <- beta;
  beta_centered[2:] <- exp(beta);
  beta_centered[1] <- p*beta_centered[2:];
  beta_aug[1] <- beta_centered[1];
  beta_centered <- Sigma*beta_centered;
}
model {
  //exponential GAM
  y ~ normal(intercept + Xs * beta_centered, sigma2);
  //P-spline prior on the differences (K-1 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  Diff*beta_aug ~ normal(0, sigma2/lambda);
}
generated quantities {
  matrix[N,K] designmat; //design matrix
  //matrix[N,K] weighted; //weighted basis functions
  real offset;
  vector[N] pred; //spline interpolant
  real edf; //effective degrees of freedom
  vector[K+1] edfvec;
  offset <- intercept;
  designmat <- Xs;
  //weighted <-  bspline(x, K, splinedegree());
  //weighted <- X .* rep_matrix(beta', rows(x));
  pred <- intercept + Xs * beta_centered;
  {
    matrix[K,K] XtX;
    XtX <- crossprod(Xs);
    edfvec[1] <- 1;
    edfvec[2:] <- diagonal(inverse_spd(XtX+lambda*crossprod(Diff)) * XtX);
    edf <- sum(edfvec);
  }
}
