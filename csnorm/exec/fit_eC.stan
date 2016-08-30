////
// Cut-site normalization model: predict mean for each provided count
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Kdiag; //number of functions in spline base for diagonal decay
  //biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  //distance bounds
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //counts : explicit
  int<lower=0> Nclose; //number of close counts modelled explicitly
  int<lower=0> counts_close[Nclose]; //value of the count
  int<lower=0> index_close[2,Nclose]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Nclose] dist_close; //genomic distance between rsites
  //
  int<lower=0> Nfar; //number of far counts modelled explicitly
  int<lower=0> counts_far[Nfar]; //value of the count
  int<lower=0> index_far[2,Nfar]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Nfar] dist_far; //genomic distance between rsites
  //
  int<lower=0> Nup; //number of upstream counts modelled explicitly
  int<lower=0> counts_up[Nup]; //value of the count
  int<lower=0> index_up[2,Nup]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Nup] dist_up; //genomic distance between rsites
  //
  int<lower=0> Ndown; //number of downstream counts modelled explicitly
  int<lower=0> counts_down[Ndown]; //value of the count
  int<lower=0> index_down[2,Ndown]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Ndown] dist_down; //genomic distance between rsites
  //
  //estimated parameters
  vector[S] log_nu; //take nu and delta directly to avoid base reconstruction
  vector[S] log_delta; 
  vector[Kdiag] beta_diag_centered; //need to build spline base
}
transformed data {
  //diag
  vector[Nclose+Nfar+Nup+Ndown] log_decay;
  vector[Nclose] log_decay_close;
  vector[Nfar] log_decay_far;
  vector[Nup] log_decay_up;
  vector[Ndown] log_decay_down;
  //means
  vector[Nclose] log_mean_cclose;
  vector[Nfar] log_mean_cfar;
  vector[Nup] log_mean_cup;
  vector[Ndown] log_mean_cdown;
  
  //decay: diagonal SCAM spline, dense
  {
    matrix[Nclose+Nfar+Nup+Ndown,Kdiag] Xdiag;
    //design matrix
    {
      vector[Nclose+Nfar+Nup+Ndown] tmpN;
      tmpN[:Nclose] = dist_close;
      tmpN[(Nclose+1):(Nclose+Nfar)] = dist_far;
      tmpN[(Nclose+Nfar+1):(Nclose+Nfar+Nup)] = dist_up;
      tmpN[(Nclose+Nfar+Nup+1):] = dist_down;
      Xdiag = bspline(log(tmpN), Kdiag, splinedegree(), log(dmin), log(dmax));
    }
    //decay
    log_decay = Xdiag * beta_diag_centered;
    log_decay_close = log_decay[:Nclose];
    log_decay_far = log_decay[(Nclose+1):(Nclose+Nfar)];
    log_decay_up = log_decay[(Nclose+Nfar+1):(Nclose+Nfar+Nup)];
    log_decay_down = log_decay[(Nclose+Nfar+Nup+1):];
  }
  //means
  log_mean_cclose = log_decay_close + (log_nu - log_delta)[index_close[1]] 
                                         + (log_nu + log_delta)[index_close[2]];
  log_mean_cfar   = log_decay_far   + (log_nu + log_delta)[index_far[1]]   
                                          + (log_nu - log_delta)[index_far[2]];
  log_mean_cup    = log_decay_up    + (log_nu + log_delta)[index_up[1]]    
                                          + (log_nu + log_delta)[index_up[2]];
  log_mean_cdown  = log_decay_down  + (log_nu - log_delta)[index_down[1]]  
                                          + (log_nu - log_delta)[index_down[2]];
}
parameters {
  real eC;  //exposure for counts
  real<lower=0> alpha;
}
model {
  //counts: Close, Far, Up, Down
  counts_close ~ neg_binomial_2_log(eC + log_mean_cclose, alpha); // Close
  counts_far   ~ neg_binomial_2_log(eC + log_mean_cfar, alpha); // Far
  counts_up    ~ neg_binomial_2_log(eC + log_mean_cup, alpha); // Up
  counts_down  ~ neg_binomial_2_log(eC + log_mean_cdown, alpha); // Down
}
generated quantities {
  vector[Nclose+Nfar+Nup+Ndown] ldec;
  ldec = log_decay;
}
