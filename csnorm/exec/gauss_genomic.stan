////
// Cut-site normalization model: normal approximation for nu and delta
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //experimental design
  int<lower=1> Dsets; //number of datasets
  int<lower=1,upper=Dsets> Biases; //number of different genomic biases to model
  int<lower=1,upper=Biases> XB[Dsets]; //XB[i]=j: dataset i has bias set j
  //genomic biases
  int<lower=4> Krow; //number of functions in spline base for row biases
  int<lower=1> SD; //number of cut sites across all datasets
  int<lower=1,upper=SD+1> bbegin[Dsets+1]; //bbegin[i]=j: dataset i starts at j
  vector[SD] cutsitesD; //cut site locations, all data
  int<lower=0> rejoined[SD];
  int<lower=0> danglingL[SD];
  int<lower=0> danglingR[SD];
  //reduced count sums
  vector[SD] eta_hat_L;
  vector[SD] eta_hat_R;
  //standard deviations
  vector<lower=0>[SD] sd_L;
  vector<lower=0>[SD] sd_R;
  //dispersion
  real<lower=0> alpha;
  //stiffnesses
  real<lower=0> lambda_nu[Biases];
  real<lower=0> lambda_delta[Biases];
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(SD)] XDrow_w;
  int XDrow_v[nnz(SD)];
  int XDrow_u[SD+1];
  row_vector[Krow] prowD[Dsets];
  //scaling factor for genomic lambdas
  real lfac;

  //Bias GAM spline, sparse
  {
    vector[SD] row_weightsD;
    int nnzs[Dsets+1];
    row_weightsD = rep_vector(1, SD);
    nnzs[1] = 1;
    for (i in 1:Dsets) nnzs[i+1] = nnzs[i]+nnz(bbegin[i+1]-bbegin[i]);
    
    #compute design matrix and projector
    for (d in 1:Dsets) {
      int S;
      S = bbegin[d+1]-bbegin[d];
      {
        vector[S] cutsites;
        vector[nnz(S)] Xrow_w;
        int Xrow_v[nnz(S)];
        int Xrow_u[S+1];
        vector[S] row_weights;
        row_vector[Krow] prow;
        cutsites = cutsitesD[bbegin[d]:(bbegin[d+1]-1)];
        row_weights = row_weightsD[bbegin[d]:(bbegin[d+1]-1)];
        #include "sparse_spline_construction.stan"
        XDrow_w[nnzs[d]:(nnzs[d+1]-1)] = Xrow_w;
        for (i in 1:size(Xrow_v)) Xrow_v[i] = Xrow_v[i] + Krow*(d-1);
        XDrow_v[nnzs[d]:(nnzs[d+1]-1)] = Xrow_v;
        prowD[d] = prow;
      }
    }
    XDrow_u[1] = 1;
    for (i in 1:SD) XDrow_u[i+1] = XDrow_u[i]+nnz(1);
  }
  
  //scaling factor
  lfac = 30000*Krow/(max(cutsitesD)-min(cutsitesD));
  
}
parameters {
  //exposures
  real eC[Dsets];
  real eRJ[Dsets];
  real eDE[Dsets];
  //spline parameters
  vector[Krow-1] beta_nu[Biases];
  vector[Krow-1] beta_delta[Biases];
}
transformed parameters {
  //nu
  vector[SD] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff[Dsets]; //2nd order difference on beta_nu_aug
  //delta
  vector[SD] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff[Dsets]; //2nd order difference on beta_delta_aug
  //means
  vector[SD] log_mean_DL;
  vector[SD] log_mean_DR;
  vector[SD] log_mean_RJ;
  //
  vector[SD] log_mean_cleft;
  vector[SD] log_mean_cright;

  //nu
  {
    vector[Dsets*Krow] beta_nu_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_nu_aug;
      vector[Krow] tmp;
      beta_nu_aug[1] = 0;
      beta_nu_aug[2:] = beta_nu[XB[d]];
      tmp = beta_nu_aug - (prowD[d] * beta_nu_aug) * rep_vector(1,Krow);
      beta_nu_centered[((d-1)*Krow+1):(d*Krow)] = tmp;
      beta_nu_diff[d] = tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_nu = csr_matrix_times_vector(SD, Dsets*Krow, XDrow_w, XDrow_v, XDrow_u, beta_nu_centered);
  }
  
  //delta
  {
    vector[Dsets*Krow] beta_delta_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_delta_aug;
      vector[Krow] tmp;
      beta_delta_aug[1] = 0;
      beta_delta_aug[2:] = beta_delta[XB[d]];
      tmp = beta_delta_aug - (prowD[d] * beta_delta_aug) * rep_vector(1,Krow);
      beta_delta_centered[((d-1)*Krow+1):(d*Krow)] = tmp;
      beta_delta_diff[d] = tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_delta = csr_matrix_times_vector(SD, Dsets*Krow, XDrow_w, XDrow_v, XDrow_u, beta_delta_centered);
  }
  
  //means
  {
    vector[SD] counts_exposure;
    //biases
    log_mean_RJ = log_nu;
    log_mean_DL = log_nu + log_delta;
    log_mean_DR = log_nu - log_delta;
    //add exposures
    for (d in 1:Dsets) {
      log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)] + eRJ[d];
      log_mean_DL[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DL[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_DR[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DR[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      counts_exposure[bbegin[d]:(bbegin[d+1]-1)] = rep_vector(eC[d],bbegin[d+1]-bbegin[d]);
    }
    //exact counts  
    log_mean_cleft  = counts_exposure + log_nu + log_delta;
    log_mean_cright = counts_exposure + log_nu - log_delta;
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts
  //grouping reduces the number of likelihoods from S-1 to G, so reweighting is
  //needed for a proper estimation of the lambdas
  eta_hat_L ~ normal(log_mean_cleft, sd_L);
  eta_hat_R ~ normal(log_mean_cright, sd_R);
  
  //// prior
  log_nu ~ cauchy(0, 1); //give high probability to [0.5:2]
  log_delta ~ cauchy(0,1);
  for (d in 1:Dsets) {
    beta_nu_diff[d] ~ normal(0,1/(lfac*lambda_nu[XB[d]]));
    beta_delta_diff[d] ~ normal(0,1/(lfac*lambda_delta[XB[d]]));
    beta_nu_diff[d] ~ double_exponential(0,10/(lfac*lambda_nu[XB[d]]));
    beta_delta_diff[d] ~ double_exponential(0,10/(lfac*lambda_delta[XB[d]]));
  }
}
