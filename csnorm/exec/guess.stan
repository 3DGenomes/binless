////
// Cut-site normalization model: initial guess
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //genomic biases
  int<lower=4> Krow; //number of functions in spline base for row biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  int<lower=0> rejoined[S];
  int<lower=0> danglingL[S];
  int<lower=0> danglingR[S];
  //count sums
  int<lower=0> counts_sum_left[S];
  int<lower=0> counts_sum_right[S];
  //stiffnesses
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;
  //rescale lambdas
  real lnu;
  real ldelta;

  row_weights = rep_vector(1, S);
  #compute design matrix and projector
  #include "sparse_spline_construction.stan"
  
  {
    real tmp;
    tmp = 30000*Krow/(max(cutsites)-min(cutsites));
    lnu = lambda_nu*tmp;
    ldelta = lambda_delta*tmp;
  }
}
parameters {
  //exposures
  real eC;
  real eRJ;
  real eDE;
  //spline parameters
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  //dispersion
  real<lower=0> alpha;
}
transformed parameters {
  //nu
  vector[S] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff; //2nd order difference on beta_nu_aug
  //delta
  vector[S] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff; //2nd order difference on beta_delta_aug
  //means
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  vector[S] log_mean_RJ;
  //
  vector[S] log_mean_cleft;
  vector[S] log_mean_cright;

  //nu
  {
    vector[Krow] beta_nu_aug;
    vector[Krow] beta_nu_centered;
    beta_nu_aug[1] = 0;
    beta_nu_aug[2:] = beta_nu;
    beta_nu_centered = beta_nu_aug - (prow * beta_nu_aug) * rep_vector(1,Krow);
    log_nu = csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
    beta_nu_diff = beta_nu_centered[:(Krow-2)]-2*beta_nu_centered[2:(Krow-1)]+beta_nu_centered[3:];
  }
  //delta
  {
    vector[Krow] beta_delta_aug;
    vector[Krow] beta_delta_centered;
    beta_delta_aug[1] = 0;
    beta_delta_aug[2:] = beta_delta;
    beta_delta_centered = beta_delta_aug - (prow * beta_delta_aug) * rep_vector(1,Krow);
    log_delta = csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
    beta_delta_diff = beta_delta_centered[:(Krow-2)]-2*beta_delta_centered[2:(Krow-1)]+beta_delta_centered[3:];
  }

  //means
  {
    //biases
    log_mean_RJ = log_nu + eRJ;
    log_mean_DL = log_nu + eDE + log_delta;
    log_mean_DR = log_nu + eDE - log_delta;
    //exact counts  
    log_mean_cleft  = eC + log_nu + log_delta;
    log_mean_cright = eC + log_nu - log_delta;
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts
  increment_log_prob((S-1) * neg_binomial_2_log_log(counts_sum_left, log_mean_cleft, alpha));
  increment_log_prob((S-1) * neg_binomial_2_log_log(counts_sum_right, log_mean_cright, alpha));
  
  //// prior
  log_nu ~ cauchy(0, 1); //give high probability to [0.5:2]
  log_delta ~ cauchy(0,1);
  beta_nu_diff ~ normal(0,1/lnu);
  beta_delta_diff ~ normal(0,1/ldelta);
  beta_nu_diff ~ double_exponential(0,10/lnu);
  beta_delta_diff ~ double_exponential(0,10/ldelta);
}
