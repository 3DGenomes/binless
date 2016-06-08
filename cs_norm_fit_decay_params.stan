////
// Cut-site normalization model: retrieve log_decay, eC and beta_diag_centered
////
functions {
  ///////////
  //BEGIN internal use
  //
  //generate [imin, imin+1, ... imax]
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
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
  
  matrix bspline(vector x, int K, int q, real xmin, real xmax) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(xmax-xmin)/(K-q); //make it slightly larger
    t <- xmin - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }

}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Kdiag; //number of functions in spline base for diagonal decay
  //distance bounds
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //input data
  int<lower=Kdiag+1> N;
  vector[N] dist;
  vector[N] fij;
  vector<lower=1>[N] ncounts;
}
transformed data {
  //diagonal SCAM spline, dense, exact
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag] pdiag;
  vector[N] diag_weights;
  
  //diagonal SCAM spline, dense, exact
  {
    Xdiag <- bspline(log(dist), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    diag_weights <- ncounts/mean(ncounts);
    pdiag <- diag_weights' * Xdiag;
    pdiag <- pdiag / (pdiag * rep_vector(1,Kdiag));
  }
}
parameters {
  //exposures
  real eC;
  positive_ordered[Kdiag-1] beta_diag;
  //SD
  real<lower=0> sigma;
  //length scales
  real<lower=0> lambda_diag;
}
transformed parameters {
  //diag
  vector[N] log_decay;
  vector[Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff;
  {
    vector[Kdiag] beta_diag_aug;
    real epsilon;
    epsilon <- -1; //decreasing spline
    beta_diag_aug[1] <- 0;
    beta_diag_aug[2:] <- beta_diag;
    beta_diag_centered <- epsilon * (beta_diag_aug - (pdiag * beta_diag_aug) * rep_vector(1, Kdiag));
    log_decay <- Xdiag * beta_diag_centered;
    beta_diag_diff <- beta_diag_centered[:(Kdiag-2)]-2*beta_diag_centered[2:(Kdiag-1)]+beta_diag_centered[3:];
  }
}
model {
  //// Exact likelihoods
  fij ~ normal(eC+log_decay, sigma);
  
  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_diag_diff ~ normal(0, 1./lambda_diag);
}
