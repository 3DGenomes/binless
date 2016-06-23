////
// Cut-site normalization model: retrieve log_decay, eC and beta_diag_centered
////
functions {
  #include "common_functions.stan"
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
  
  #include "scam_spline_construction.stan"
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
