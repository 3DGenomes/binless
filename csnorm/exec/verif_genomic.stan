////
// Cut-site normalization model: simplified fit for nu and delta
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
  //counts : explicit
  int<lower=0> N; //number of counts
  int<lower=1,upper=S> cidx[2,N]; //indices of rsite pairs
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
  //fitted parameters
  vector[S] log_nu_init;
  vector[S] log_delta_init;
  vector[N] log_decay;
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;
  //scaling factor for genomic lambdas
  real lfac;
  
  row_weights <- rep_vector(1, S);
  #compute design matrix and projector
  #include "sparse_spline_construction.stan"
  
  //scaling factor
  lfac <- 30000*Krow/(max(cutsites)-min(cutsites));
  
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
  //stiffnesses
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
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
  vector[N] log_mean_cclose1;
  vector[N] log_mean_cfar1;
  vector[N] log_mean_cup1;
  vector[N] log_mean_cdown1;
  //
  vector[N] log_mean_cclose2;
  vector[N] log_mean_cfar2;
  vector[N] log_mean_cup2;
  vector[N] log_mean_cdown2;
  
  //nu
  {
    vector[Krow] beta_nu_aug;
    vector[Krow] beta_nu_centered;
    beta_nu_aug[1] <- 0;
    beta_nu_aug[2:] <- beta_nu;
    beta_nu_centered <- beta_nu_aug - (prow * beta_nu_aug) * rep_vector(1,Krow);
    log_nu <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
    beta_nu_diff <- beta_nu_centered[:(Krow-2)]-2*beta_nu_centered[2:(Krow-1)]+beta_nu_centered[3:];
  }
  //delta
  {
    vector[Krow] beta_delta_aug;
    vector[Krow] beta_delta_centered;
    beta_delta_aug[1] <- 0;
    beta_delta_aug[2:] <- beta_delta;
    beta_delta_centered <- beta_delta_aug - (prow * beta_delta_aug) * rep_vector(1,Krow);
    log_delta <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
    beta_delta_diff <- beta_delta_centered[:(Krow-2)]-2*beta_delta_centered[2:(Krow-1)]+beta_delta_centered[3:];
  }

  //means
  {
    //biases
    log_mean_RJ <- log_nu + eRJ;
    log_mean_DL <- log_nu + eDE + log_delta;
    log_mean_DR <- log_nu + eDE - log_delta;
    //exact counts  
    log_mean_cclose1 <- eC + log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu_init + log_delta_init)[cidx[2]];
    log_mean_cfar1   <- eC + log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu_init - log_delta_init)[cidx[2]];
    log_mean_cup1    <- eC + log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu_init + log_delta_init)[cidx[2]];
    log_mean_cdown1  <- eC + log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu_init - log_delta_init)[cidx[2]];
    //
    log_mean_cclose2 <- eC + log_decay + (log_nu_init - log_delta_init)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cfar2   <- eC + log_decay + (log_nu_init + log_delta_init)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
    log_mean_cup2    <- eC + log_decay + (log_nu_init + log_delta_init)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cdown2  <- eC + log_decay + (log_nu_init - log_delta_init)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts: Close, Far, Up, Down
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_close, log_mean_cclose1, alpha));
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_far, log_mean_cfar1, alpha));
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_up, log_mean_cup1, alpha));
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_down, log_mean_cdown1, alpha));
  //
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_close, log_mean_cclose2, alpha));
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_far, log_mean_cfar2, alpha));
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_up, log_mean_cup2, alpha));
  increment_log_prob(0.5*neg_binomial_2_log_log(counts_down, log_mean_cdown2, alpha));
  
  //// prior
  log_nu ~ cauchy(0, 1); //give high probability to [0.5:2]
  log_delta ~ cauchy(0,1);
  beta_nu_diff ~ normal(0,1/(lfac*lambda_nu));
  beta_delta_diff ~ normal(0,1/(lfac*lambda_delta));
  beta_nu_diff ~ double_exponential(0,10/(lfac*lambda_nu));
  beta_delta_diff ~ double_exponential(0,10/(lfac*lambda_delta));
}
