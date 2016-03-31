////
// GAM : y ~ s(X) using cubic splines with log-link
// The design matrix is the original spline base
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
  matrix design_spline(vector x, int K, int q) {
    return bspline(x, K, q);
  }
  
  // difference operator for that design matrix
  matrix difference_op(int K, int d) {
    return difference_matrix_sqrt(K,d);
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
  matrix[N,K] Xs; //design matrix for spline
  matrix[K-difforder(),K] Diff; //difference operator for penalty
  vector[K] p; //projector on centering constraint
  Xs <- design_spline(x, K, splinedegree());
  Diff <- difference_op(K, difforder());
  p <- Xs' * rep_vector(1,N);
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
  y ~ neg_binomial_2(exp(intercept + Xs * beta_centered), alpha);
  //P-spline prior on the differences (K-1 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  Diff*beta_centered ~ normal(0, 1./(alpha*lambda));
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
  pred <- exp(intercept + Xs * beta_centered);
  {
    matrix[K,K] XtX;
    XtX <- crossprod(Xs);
    edfvec[1] <- 1;
    edfvec[2:] <- diagonal(inverse_spd(XtX+lambda*crossprod(Diff)) * XtX);
    edf <- sum(edfvec);
  }
}
