////
// Cut-site normalization model: retrieve eRJ eDE log_nu and log_delta
//                               at constant stiffness
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int<lower=0> Krow; //number of functions in spline base for row biases
  //biases
  int<lower=1> S;
  vector[S] cutsites;
  vector[S] log_mean_RJ;
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
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
  
  row_weights <- rep_vector(1, S);
  #compute design matrix and projector
  #include "sparse_spline_construction.stan"

  {
    real tmp;
    tmp <- 30000*Krow/(max(cutsites)-min(cutsites));
    lnu <- lambda_nu*tmp;
    ldelta <- lambda_delta*tmp;
  }
}
parameters {
  //exposures
  real eRJ;
  real eDE;
  //spline coefs
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  //SD
  real<lower=0,upper=1> sigma;
}
transformed parameters {
  //nu
  vector[S] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff; //2nd order difference on beta_nu_aug
  //delta
  vector[S] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff; //2nd order difference on beta_delta_aug
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
}
model {
  //Likelihoods
  log_mean_RJ ~ normal(eRJ + log_nu, sigma);
  log_mean_DL ~ normal(eDE  + log_nu + log_delta, sigma);
  log_mean_DR ~ normal(eDE  + log_nu - log_delta, sigma);
  //P-spline prior on the 2nd order differences (Krow-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_nu_diff ~ normal(0, 1/lnu);
  beta_delta_diff ~ normal(0, 1/ldelta);
}
