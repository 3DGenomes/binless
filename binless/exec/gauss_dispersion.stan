////
// Cut-site normalization model: fit dispersion and exposures only
////
////////////////////////////////////////////////////////////////
data {
  //experimental design
  int<lower=1> Dsets; //number of datasets
  int<lower=1,upper=Dsets> Biases; //number of different genomic biases to model
  int<lower=1,upper=Dsets> Decays; // number of different diagonal decays to model
  int<lower=1,upper=Biases> XB[Dsets]; //XB[i]=j: dataset i has bias set j
  int<lower=1,upper=Decays> XD[Dsets]; //XD[i]=j: dataset i has decay set j
  //genomic biases
  int<lower=1> SD; //number of cut sites across all datasets
  int<lower=1,upper=SD+1> bbegin[Dsets+1]; //bbegin[i]=j: dataset i starts at j
  vector[SD] cutsitesD; //cut site locations, all data
  int rejoined[SD];
  int danglingL[SD];
  int danglingR[SD];
  //counts : explicit
  int<lower=0> N; //number of counts
  int<lower=1,upper=N+1> cbegin[Dsets+1]; //cbegin[i]=j: dataset i starts at j
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
  real<lower=0> weight[Dsets]; //in case of subsampling
  //parameters: spline and means
  vector[SD] log_iota;
  vector[SD] log_rho;
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
}
parameters {
  //exposures
  real eC_sup[Dsets];
  real eRJ[Dsets];
  real eDE[Dsets];
  //dispersion
  real<lower=0> alpha;
}
transformed parameters {
  //means
  vector[SD] log_mean_DL;
  vector[SD] log_mean_DR;
  vector[SD] log_mean_RJ;
    
  //means
  {
    //biases
    log_mean_RJ = (log_iota + log_rho)/2;
    log_mean_DL = log_iota;
    log_mean_DR = log_rho;
    //add exposures
    for (d in 1:Dsets) {
      log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)] + eRJ[d];
      log_mean_DL[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DL[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_DR[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DR[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
    }
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts: Close, Far, Up, Down
  for (d in 1:Dsets) {
    target += weight[d]*neg_binomial_2_log_lpmf(counts_close[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] + eC_sup[d], alpha);
    target += weight[d]*neg_binomial_2_log_lpmf(counts_far[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)] + eC_sup[d], alpha);
    target += weight[d]*neg_binomial_2_log_lpmf(counts_up[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cup[cbegin[d]:(cbegin[d+1]-1)] + eC_sup[d], alpha);
    target += weight[d]*neg_binomial_2_log_lpmf(counts_down[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)] + eC_sup[d], alpha);
  }
}
