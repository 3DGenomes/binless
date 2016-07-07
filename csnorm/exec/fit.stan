////
// Cut-site normalization model: fit available data, partial mean-field approximation
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
  int rejoined[S];
  int danglingL[S];
  int danglingR[S];
  //decay bias
  int<lower=4> Kdiag; //number of functions in spline base for diagonal decay
  real<lower=0> dmin; //distance bounds for spline
  real<lower=dmin> dmax;
  //counts : explicit
  int<lower=0> N; //number of counts
  int<lower=1,upper=S> cidx[2,N]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[N] dist; //genomic distance between rsites
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;
  //diagonal SCAM spline, dense
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag] pdiag;
  
  row_weights <- rep_vector(3, S); //dangling L/R + rejoined
  for (i in 1:N) {
    row_weights[cidx[1,i]] <- row_weights[cidx[1,i]] + 1;
    row_weights[cidx[2,i]] <- row_weights[cidx[2,i]] + 1;
  }
  #compute design matrix and projector
  #include "sparse_spline_construction.stan"

  //diagonal SCAM spline, dense, exact and mean field model
  #include "scam_spline_construction.stan"
}
parameters {
  //exposures
  real eC;
  real eRJ;
  real eDE;
  //spline parameters
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  positive_ordered[Kdiag-1] beta_diag;
  //dispersion
  real<lower=0> alpha;
  //length scales
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
  real<lower=0> lambda_diag;
}
transformed parameters {
  //nu
  vector[S] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff; //2nd order difference on beta_nu_aug
  //delta
  vector[S] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff; //2nd order difference on beta_delta_aug
  //diag
  vector[N] log_decay;
  vector[Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff;
  //means
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  vector[S] log_mean_RJ;
  //
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  
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
  //decay
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

  //means
  {
    //biases
    log_mean_RJ <- log_nu + eRJ;
    log_mean_DL <- log_nu + eDE + log_delta;
    log_mean_DR <- log_nu + eDE - log_delta;
    //exact counts  
    log_mean_cclose <- eC + log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cfar   <- eC + log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
    log_mean_cup    <- eC + log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cdown  <- eC + log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts: Close, Far, Up, Down
  counts_close ~ neg_binomial_2_log(log_mean_cclose, alpha); // Close
  counts_far   ~ neg_binomial_2_log(log_mean_cfar, alpha); // Far
  counts_up    ~ neg_binomial_2_log(log_mean_cup, alpha); // Up
  counts_down  ~ neg_binomial_2_log(log_mean_cdown, alpha); // Down
  
  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_nu_diff ~ normal(0, lambda_nu);
  beta_delta_diff ~ normal(0, lambda_delta);
  beta_diag_diff ~ normal(0, lambda_diag);
  //weak lasso on bias splines (non-bayesian)
  beta_nu_diff ~ double_exponential(0,10*lambda_nu);
  beta_delta_diff ~ double_exponential(0,10*lambda_delta);

  //cauchy hyperpriors for genomic biases
  log_nu ~ cauchy(0,1);
  log_delta ~ cauchy(0,1);

  //half cauchy hyperpriors for lengthscales
  lambda_nu ~ cauchy(0,0.1);
  lambda_delta ~ cauchy(0,0.1);
  lambda_diag ~ cauchy(0,0.1);
}
generated quantities {
  real deviance;
  real deviance_null;
  real deviance_proportion_explained;
  #
  real rejoined_deviance;
  real rejoined_deviance_null;
  real rejoined_deviance_proportion_explained;
  #
  real dangling_deviance;
  real dangling_deviance_null;
  real dangling_deviance_proportion_explained;
  #
  real count_deviance;
  real count_deviance_null;
  real count_deviance_proportion_explained;
  #deviances
  count_deviance <- neg_binomial_2_log_deviance(counts_close, log_mean_cclose, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_far, log_mean_cfar, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_up, log_mean_cup, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_down, log_mean_cdown, alpha, rep_vector(1,1));
  rejoined_deviance <- neg_binomial_2_log_deviance(rejoined, log_mean_RJ, alpha, rep_vector(1,1));
  dangling_deviance <- neg_binomial_2_log_deviance(danglingL, log_mean_DL, alpha, rep_vector(1,1)) +
                       neg_binomial_2_log_deviance(danglingR, log_mean_DR, alpha, rep_vector(1,1));
  deviance <- rejoined_deviance + dangling_deviance + count_deviance;
  #null deviances
  rejoined_deviance_null <- neg_binomial_2_log_deviance(rejoined, rep_vector(eRJ,1), alpha, rep_vector(1,1));
  dangling_deviance_null <- neg_binomial_2_log_deviance(danglingL, rep_vector(eDE,1), alpha, rep_vector(1,1)) +
                            neg_binomial_2_log_deviance(danglingR, rep_vector(eDE,1), alpha, rep_vector(1,1));
  count_deviance_null <- neg_binomial_2_log_deviance(counts_close, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_far, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_up, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_down, rep_vector(eC,1), alpha, rep_vector(1,1));
  deviance_null <- rejoined_deviance_null + dangling_deviance_null + count_deviance_null;
  #proportions explained
  rejoined_deviance_proportion_explained <- 100*(rejoined_deviance_null - rejoined_deviance)/rejoined_deviance_null;
  dangling_deviance_proportion_explained <- 100*(dangling_deviance_null - dangling_deviance)/dangling_deviance_null;
  count_deviance_proportion_explained <- 100*(count_deviance_null - count_deviance)/count_deviance_null;
  deviance_proportion_explained <- 100*(deviance_null - deviance)/deviance_null;
}
