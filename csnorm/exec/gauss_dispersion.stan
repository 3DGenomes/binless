////
// Cut-site normalization model: fit dispersion and exposures only
////
functions {
  #include "common_functions.stan"
}
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
  //decay bias
  int<lower=4> Kdiag; //number of functions in spline base for diagonal decay
  real<lower=0> dmin; //distance bounds for spline
  real<lower=dmin> dmax;
  //counts : explicit
  int<lower=0> N; //number of counts
  int<lower=1,upper=N+1> cbegin[Dsets+1]; //cbegin[i]=j: dataset i starts at j
  int<lower=1,upper=SD> cidx[2,N]; //indices of rsite pairs, referring to vector[SD] cutsitesD
  vector<lower=dmin,upper=dmax>[N] dist; //genomic distance between rsites
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
  real<lower=0> weight[Dsets]; //in case of subsampling
  //parameters: spline
  vector[SD] log_iota;
  vector[SD] log_rho;
  positive_ordered[Kdiag-1] beta_diag[Decays];
}
transformed data {
  //diagonal SCAM spline, dense
  matrix[N,Dsets*Kdiag] Xdiag;
  row_vector[Kdiag] pdiagD[Dsets];
  //diag
  vector[N] ldec;
  vector[Dsets*Kdiag] beta_diag_centered;
  
  //diagonal SCAM spline, dense
  {
    {
      matrix[N,Kdiag] tmpXdiag;
      tmpXdiag = bspline(log(dist), Kdiag, splinedegree(), log(dmin), log(dmax));
      Xdiag = rep_matrix(0, N, Dsets*Kdiag);
      for (d in 1:Dsets) {
        Xdiag[cbegin[d]:(cbegin[d+1]-1),((d-1)*Kdiag+1):(d*Kdiag)] = tmpXdiag[cbegin[d]:(cbegin[d+1]-1),:];
      }
    }
      
    //projector for diagonal (SCAM)
    for (d in 1:Dsets) {
      int sz;
      sz = cbegin[d+1]-cbegin[d];
      {
        row_vector[sz] diag_weights;
        row_vector[Kdiag] pdiag;
        diag_weights = rep_row_vector(1,sz);
        diag_weights = diag_weights/mean(diag_weights);
        pdiag = diag_weights * Xdiag[cbegin[d]:(cbegin[d+1]-1),((d-1)*Kdiag+1):(d*Kdiag)];
        pdiagD[d] = pdiag / (pdiag * rep_vector(1,Kdiag));
      }
    }
  }
  
  //decay
  for (d in 1:Dsets) {
    vector[Kdiag] beta_diag_aug;
    vector[Kdiag] tmp;
    real epsilon;
    epsilon = -1; //decreasing spline
    beta_diag_aug[1] = 0;
    beta_diag_aug[2:] = beta_diag[XD[d]];
    tmp = epsilon * (beta_diag_aug - (pdiagD[d] * beta_diag_aug) * rep_vector(1, Kdiag));
    beta_diag_centered[((d-1)*Kdiag+1):(d*Kdiag)] = tmp;
  }
  ldec = Xdiag * beta_diag_centered;

}
parameters {
  //exposures
  real eC[Dsets];
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
  //
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
    
  //means
  {
    //biases
    log_mean_RJ = (log_iota + log_rho)/2;
    log_mean_DL = log_iota;
    log_mean_DR = log_rho;
    //exact counts  
    log_mean_cclose = ldec + log_rho[cidx[1]]  + log_iota[cidx[2]];
    log_mean_cfar   = ldec + log_iota[cidx[1]] + log_rho[cidx[2]];
    log_mean_cup    = ldec + log_iota[cidx[1]] + log_iota[cidx[2]];
    log_mean_cdown  = ldec + log_rho[cidx[1]]  + log_rho[cidx[2]];
    //add exposures
    for (d in 1:Dsets) {
      log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)] + eRJ[d];
      log_mean_DL[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DL[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_DR[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DR[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] = log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
      log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)]   = log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
      log_mean_cup[cbegin[d]:(cbegin[d+1]-1)]    = log_mean_cup[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
      log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)]  = log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
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
    target += weight[d]*neg_binomial_2_log_lpmf(counts_close[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)], alpha);
    target += weight[d]*neg_binomial_2_log_lpmf(counts_far[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)], alpha);
    target += weight[d]*neg_binomial_2_log_lpmf(counts_up[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cup[cbegin[d]:(cbegin[d+1]-1)], alpha);
    target += weight[d]*neg_binomial_2_log_lpmf(counts_down[cbegin[d]:(cbegin[d+1]-1)] | log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)], alpha);
  }
}
generated quantities {
  vector[N] log_decay;
  log_decay = ldec;
}

