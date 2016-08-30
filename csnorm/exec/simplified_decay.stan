////
// Cut-site normalization model: simplified fit for decay
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //genomic biases
  int<lower=1> S; //number of cut sites
  int rejoined[S];
  int danglingL[S];
  int danglingR[S];
  //estimated nu and delta
  vector[S] log_nu;
  vector[S] log_delta;
  //decay bias
  int<lower=4> Kdiag; //number of functions in spline base for diagonal decay
  real<lower=0> dmin; //distance bounds for spline
  real<lower=dmin> dmax;
  //count sums
  int<lower=1> N;
  int<lower=0> counts_sum[N];
  vector<lower=0>[N] weight;
  vector<lower=0>[N] dist;
  //genomic bias sums
  vector[N] log_genomic_sum;
}
transformed data {
  //diagonal SCAM spline, dense
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag] pdiag;
  row_vector[N] diag_weights;
  diag_weights = weight';
  
  //diagonal SCAM spline, dense, exact and mean field model
  #include "scam_spline_construction.stan"
}
parameters {
  //exposures
  real eC;
  real eRJ;
  real eDE;
  //spline parameters
  positive_ordered[Kdiag-1] beta_diag;
  //dispersion
  real<lower=0> alpha;
  //length scales
  real<lower=0> lambda_diag;
}
transformed parameters {
  //diag
  vector[N] log_decay;
  vector[Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff;
  //means
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  vector[S] log_mean_RJ;
  //
  vector[N] log_mean_counts;
  
  //decay
  {
    vector[Kdiag] beta_diag_aug;
    real epsilon;
    epsilon = -1; //decreasing spline
    beta_diag_aug[1] = 0;
    beta_diag_aug[2:] = beta_diag;
    beta_diag_centered = epsilon * (beta_diag_aug - (pdiag * beta_diag_aug) * rep_vector(1, Kdiag));
    log_decay = Xdiag * beta_diag_centered;
    beta_diag_diff = beta_diag_centered[:(Kdiag-2)]-2*beta_diag_centered[2:(Kdiag-1)]+beta_diag_centered[3:];
  }

  //means
  {
    //biases
    log_mean_RJ = log_nu + eRJ;
    log_mean_DL = log_nu + eDE + log_delta;
    log_mean_DR = log_nu + eDE - log_delta;
    //exact counts  
    log_mean_counts  = eC + log_decay + log_genomic_sum;
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts
  for (i in 1:N) increment_log_prob(weight[i] * neg_binomial_2_log_log(counts_sum[i], log_mean_counts[i], alpha));
  
  //// prior
  beta_diag_diff ~ normal(0,1/lambda_diag);
}
