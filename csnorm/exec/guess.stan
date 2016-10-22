////
// Cut-site normalization model: initial guess
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
  int rejoined[SD];
  int danglingL[SD];
  int danglingR[SD];
  //count sums
  int<lower=0> counts_sum_left[SD];
  int<lower=0> counts_sum_right[SD];
  //stiffnesses
  real<lower=0> lambda_iota;
  real<lower=0> lambda_rho;
}
transformed data {
  //bias spline, sparse (iota and rho have the same design)
  vector[nnz(SD)] XDrow_w;
  int XDrow_v[nnz(SD)];
  int XDrow_u[SD+1];
  row_vector[Krow] prowD[Dsets];
  //rescale lambdas
  real liota;
  real lrho;

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
  
  {
    real tmp;
    tmp = 30000*Krow/(max(cutsitesD)-min(cutsitesD));
    liota = lambda_iota*tmp;
    lrho = lambda_rho*tmp;
  }
}
parameters {
  //exposures
  real eC[Dsets];
  real eRJ[Dsets];
  real eDE[Dsets];
  //spline parameters
  vector[Krow-1] beta_iota[Biases];
  vector[Krow-1] beta_rho[Biases];
  //dispersion
  real<lower=0> alpha;
}
transformed parameters {
  //iota
  vector[SD] log_iota; // log(iota)
  vector[Krow-2] beta_iota_diff[Dsets]; //2nd order difference on beta_iota_aug
  //rho
  vector[SD] log_rho; // log(rho)
  vector[Krow-2] beta_rho_diff[Dsets]; //2nd order difference on beta_rho_aug
  //means
  vector[SD] log_mean_DL;
  vector[SD] log_mean_DR;
  vector[SD] log_mean_RJ;
  //
  vector[SD] log_mean_cleft;
  vector[SD] log_mean_cright;

  //iota
  {
    vector[Dsets*Krow] beta_iota_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_iota_aug;
      vector[Krow] tmp;
      beta_iota_aug[1] = 0;
      beta_iota_aug[2:] = beta_iota[XB[d]];
      tmp = beta_iota_aug - (prowD[d] * beta_iota_aug) * rep_vector(1,Krow);
      beta_iota_centered[((d-1)*Krow+1):(d*Krow)] = tmp;
      beta_iota_diff[d] = tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_iota = csr_matrix_times_vector(SD, Dsets*Krow, XDrow_w, XDrow_v, XDrow_u, beta_iota_centered);
  }
  
  //rho
  {
    vector[Dsets*Krow] beta_rho_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_rho_aug;
      vector[Krow] tmp;
      beta_rho_aug[1] = 0;
      beta_rho_aug[2:] = beta_rho[XB[d]];
      tmp = beta_rho_aug - (prowD[d] * beta_rho_aug) * rep_vector(1,Krow);
      beta_rho_centered[((d-1)*Krow+1):(d*Krow)] = tmp;
      beta_rho_diff[d] = tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_rho = csr_matrix_times_vector(SD, Dsets*Krow, XDrow_w, XDrow_v, XDrow_u, beta_rho_centered);
  }
  
  //means
  {
    //biases
    log_mean_RJ = (log_iota + log_rho)/2;
    log_mean_DL = log_iota;
    log_mean_DR = log_rho;
    //exact counts  
    log_mean_cleft  = log_iota;
    log_mean_cright = log_rho;
    //add exposures
    for (d in 1:Dsets) {
      log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)] + eRJ[d];
      log_mean_DL[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DL[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_DR[bbegin[d]:(bbegin[d+1]-1)]     = log_mean_DR[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_cleft[bbegin[d]:(bbegin[d+1]-1)] = log_mean_cleft[bbegin[d]:(bbegin[d+1]-1)] + eC[d];
      log_mean_cright[bbegin[d]:(bbegin[d+1]-1)]   = log_mean_cright[bbegin[d]:(bbegin[d+1]-1)] + eC[d];
    }
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts
  for (d in 1:Dsets) {
    target += (bbegin[d+1]-bbegin[d]-1) * neg_binomial_2_log_lpmf(
      counts_sum_left[bbegin[d]:(bbegin[d+1]-1)] | log_mean_cleft[bbegin[d]:(bbegin[d+1]-1)], alpha);
    target += (bbegin[d+1]-bbegin[d]-1) * neg_binomial_2_log_lpmf(
      counts_sum_right[bbegin[d]:(bbegin[d+1]-1)] | log_mean_cright[bbegin[d]:(bbegin[d+1]-1)], alpha);
  }
  
  //// prior
  log_iota ~ cauchy(0, 1); //give high probability to [0.5:2]
  log_rho ~ cauchy(0,1);
  for (d in 1:Dsets) {
    beta_iota_diff[d] ~ normal(0,1/liota);
    beta_rho_diff[d] ~ normal(0,1/lrho);
    beta_iota_diff[d] ~ double_exponential(0,10/liota);
    beta_rho_diff[d] ~ double_exponential(0,10/lrho);
  }
}
