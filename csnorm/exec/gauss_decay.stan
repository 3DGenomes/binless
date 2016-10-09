////
// Cut-site normalization model: simplified fit for decay
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //experimental design
  int<lower=1> Dsets; //number of datasets
  int<lower=1,upper=Dsets> Decays; // number of different diagonal decays to model
  int<lower=1,upper=Decays> XD[Dsets]; //XD[i]=j: dataset i has decay set j
  //decay bias
  int<lower=4> Kdiag; //number of functions in spline base for diagonal decay
  real<lower=0> dmin; //distance bounds for spline
  real<lower=dmin> dmax;
  //count sums
  int<lower=1> N;
  int<lower=1,upper=N+1> cbegin[Dsets+1]; //cbegin[i]=j: dataset i starts at j
  vector[N] kappa_hat;
  vector<lower=0>[N] sdl;
  vector<lower=0>[N] dist;
  //dispersion
  real<lower=0> alpha;
  //length scales
  real<lower=0> lambda_diag[Decays];
}
transformed data {
  //diagonal SCAM spline, dense
  matrix[N,Dsets*Kdiag] Xdiag;
  row_vector[Kdiag] pdiagD[Dsets];
  
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
        row_vector[sz] dw;
        row_vector[Kdiag] pdiag;
        dw = rep_row_vector(1, sz);
        dw = dw/mean(dw);
        pdiag = dw * Xdiag[cbegin[d]:(cbegin[d+1]-1),((d-1)*Kdiag+1):(d*Kdiag)];
        pdiagD[d] = pdiag / (pdiag * rep_vector(1,Kdiag));
      }
    }
  }
  
}
parameters {
  //exposures
  real eC[Dsets];
  //spline parameters
  positive_ordered[Kdiag-1] beta_diag[Decays];
  //normal sd nuisance
  real<lower=0> sigma;
}
transformed parameters {
  //diag
  vector[N] log_decay;
  vector[Dsets*Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff[Dsets];
  //means
  vector[N] log_mean_counts;
  
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
    beta_diag_diff[d] = tmp[:(Kdiag-2)]-2*tmp[2:(Kdiag-1)]+tmp[3:];
  }
  log_decay = Xdiag * beta_diag_centered;
  
  //means
  {
    //exact counts  
    for (d in 1:Dsets) log_mean_counts[cbegin[d]:(cbegin[d+1]-1)] = rep_vector(eC[d],cbegin[d+1]-cbegin[d]);
    log_mean_counts  = log_decay + log_mean_counts;
  }
}
model {
  //// likelihoods
  //counts
  kappa_hat ~ normal(log_mean_counts, sdl*sigma);
  
  //// prior
  for (d in 1:Dsets) beta_diag_diff[d] ~ normal(0,1/lambda_diag[XD[d]]);
}
