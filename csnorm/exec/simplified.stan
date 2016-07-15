////
// Cut-site normalization model: simplified model with grouping of counts
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
  //count sums
  int<lower=1> G; //number of count batches
  int<lower=0> counts_sum_left[S,G];
  int<lower=0> counts_sum_right[S,G];
  //decay sums
  matrix[S,G] log_decay_sum;
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
  //scaling factor for genomic lambdas
  real lfac;
  //vectorized counts
  int csl[S*G];
  int csr[S*G];

  row_weights <- rep_vector(1, S);
  #compute design matrix and projector
  #include "sparse_spline_construction.stan"

  //diagonal SCAM spline, dense, exact and mean field model
  #include "scam_spline_construction.stan"

  //scaling factor for genomic lambdas
  lfac <- 30000*Krow/(max(cutsites)-min(cutsites));
  
  //vectorized counts
  //need to transpose to be consistent with to_vector(matrix)
  {
    int begin;
    begin<-1;
    for (i in 1:G) {
      csl[begin:(begin+S-1)] <- counts_sum_left[,i];
      csr[begin:(begin+S-1)] <- counts_sum_right[,i];
      begin<-begin+S;
    }
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
  //
  vector[S*G] log_mean_cleft;
  vector[S*G] log_mean_cright;

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
    //summed counts  
    log_mean_cleft  <- to_vector(log_decay_sum + rep_matrix(eC + log_nu + log_delta, G));
    log_mean_cright <- to_vector(log_decay_sum + rep_matrix(eC + log_nu - log_delta, G));
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts: Close, Far, Up, Down
  increment_log_prob(neg_binomial_2_log_log(counts_close, log_mean_cclose, alpha));
  increment_log_prob(neg_binomial_2_log_log(counts_far, log_mean_cfar, alpha));
  increment_log_prob(neg_binomial_2_log_log(counts_up, log_mean_cup, alpha));
  increment_log_prob(neg_binomial_2_log_log(counts_down, log_mean_cdown, alpha));
  
  //summed counts
  //grouping reduces the number of likelihoods from S-1 to G, so reweighting is
  //needed for a proper estimation of the lambdas
  increment_log_prob((S-1)/G * neg_binomial_2_log_log(csl, log_mean_cleft, alpha));
  increment_log_prob((S-1)/G * neg_binomial_2_log_log(csr, log_mean_cright, alpha));
  
  //// prior
  log_nu ~ cauchy(0, 1); //give high probability to [0.5:2]
  log_delta ~ cauchy(0,1);
  beta_nu_diff ~ normal(0,1/(lfac*lambda_nu));
  beta_delta_diff ~ normal(0,1/(lfac*lambda_delta));
  beta_nu_diff ~ double_exponential(0,10/(lfac*lambda_nu));
  beta_delta_diff ~ double_exponential(0,10/(lfac*lambda_delta));
}
