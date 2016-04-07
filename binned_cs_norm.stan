////
// GAM : N ~ s(log(distance)) + log(rejoined.1) + log(rejoined.2)
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

  // design matrix for count ~ s(log(dist)) (centered q-spline with K basis functions)
  // X is N x K-1
  matrix design_spline(vector dist, int K, int q) {
    matrix[rows(dist),K] X;
    X <- bspline(log(dist), K-q, q);
    return center(X,X);
  }
  
  // difference operator for centered spline design matrix
  matrix difference_op(vector dist, int K, int q, int d) {
    matrix[K-d,K-1] diff;
    diff <- center(bspline(log(dist), K-q, q), difference_matrix_sqrt(K,d));
    return diff;
  }
  
  //design matrix for count_k ~ sum_i b_ik (different bias coefficients for each of 4 count types) 
  matrix design_biases(matrix biases, int[] type) {
    int b;
    matrix[rows(biases), 4*cols(biases)] X;
    b <- cols(biases);
    X <- rep_matrix(0,rows(X),cols(X));
    for (i in 1:rows(X)) {
      int j;
      j <- b*(type[i]-1)+1;
      X[i,j:(j+b-1)] <- log(biases[i]);
    }
    return X;
  }
  
}
////////////////////////////////////////////////////////////////
data {
  int<lower=1> N; //number of data points
  int<lower=1> K; //number of parameters
  int<lower=1> B; //number of biases
  int count[N]; //dependent variable
  int<lower=1,upper=4> type[N]; //count type (Up Down Close Far)
  vector[N] dist; //genomic distance for this count
  matrix[N,2*B] biases; //2 biases per count
}
transformed data {
  matrix[N,K-1] X_spline; //design matrix
  matrix[N,4*2*B] X_biases;
  matrix[K-difforder(),K-1] P;
  X_spline <- design_spline(dist, K, splinedegree());
  P <- difference_op(dist, K, splinedegree(), difforder());
  X_biases <- design_biases(biases, type);
}
parameters {
  real offset;
  vector[K-1] beta_spline;
  vector[4*2*B] beta_biases;
  real<lower=0> alpha;
  real<lower=0> lambda;
}
/*transformed parameters {
  vector[N] mu;
  mu <- ;
}*/
model {
  //negative binomial GAM with log link
  count ~ neg_binomial_2_log(rep_vector(offset, N) + X_spline*beta_spline + X_biases*beta_biases, alpha);
  //P-spline prior on the differences (K-1 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  P*beta_spline ~ normal(0, 1./(alpha*lambda));
}
generated quantities {
  vector[N] pred; //mean predictive
  vector[N] decay; //spline interpolant
  vector[K+4*2*B] edfvec;
  matrix[N,4*2*B] X_b;
  real edf; //effective degrees of freedom
  pred <- exp(rep_vector(offset, N) + X_spline*beta_spline + X_biases*beta_biases);
  X_b <- X_biases;
  decay <- exp(X_spline*beta_spline);
  {
    matrix[K-1,K-1] XtX;
    XtX <- crossprod(X_spline);
    edfvec[1] <- 1;
    edfvec[2:K] <- diagonal(inverse_spd(XtX+lambda*crossprod(P)) * XtX);
    edfvec[(K+1):] <- rep_vector(1,4*2*B);
    edf <- sum(edfvec);
  }
}
