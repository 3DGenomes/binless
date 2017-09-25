////
// Cut-site normalization model: predict mean for each provided count
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int<lower=1> Dsets; //number of datasets
  int<lower=1,upper=Dsets> Decays; // number of different diagonal decays to model
  int<lower=1,upper=Decays> XD[Dsets]; //XD[i]=j: dataset i has decay set j
  int Kdiag; //number of functions in spline base for diagonal decay
  //biases
  int<lower=1> SD; //number of cut sites
  vector[SD] cutsitesD; //cut site locations
  //distance bounds
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //counts
  int<lower=0> N; //number of counts modelled explicitly
  int<lower=1,upper=N+1> cbegin[Dsets+1]; //cbegin[i]=j: dataset i starts at j
  int<lower=1,upper=SD> cidx[2,N]; //indices of rsite pairs, referring to vector[SD] cutsitesD
  vector<lower=dmin,upper=dmax>[N] dist; //genomic distance between rsites
  //estimated parameters
  real eC[Dsets];  //exposure for counts
  vector[SD] log_iota; //take iota and rho directly to avoid base reconstruction
  vector[SD] log_rho; 
  vector[Dsets*Kdiag] beta_diag_centered; //need to build spline base
}
parameters {}
model {}
generated quantities {
  //diag
  vector[N] log_decay;
  //means
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  
  //decay: diagonal SCAM spline, dense
  {
    matrix[N,Dsets*Kdiag] Xdiag;
    //design matrix
    {
      matrix[N,Kdiag] tmpXdiag;
      tmpXdiag = bspline(dist, Kdiag, splinedegree(), dmin, dmax);
      Xdiag = rep_matrix(0, N, Dsets*Kdiag);
      for (d in 1:Dsets) {
        Xdiag[cbegin[d]:(cbegin[d+1]-1),((d-1)*Kdiag+1):(d*Kdiag)] = tmpXdiag[cbegin[d]:(cbegin[d+1]-1),:];
      }
    }
    //decay
    log_decay = Xdiag * beta_diag_centered;
  }
  
  //means
  log_mean_cclose = log_decay + log_rho[cidx[1]]  + log_iota[cidx[2]];
  log_mean_cfar   = log_decay + log_iota[cidx[1]] + log_rho[cidx[2]];
  log_mean_cup    = log_decay + log_iota[cidx[1]] + log_iota[cidx[2]];
  log_mean_cdown  = log_decay + log_rho[cidx[1]]  + log_rho[cidx[2]];
  //add exposures
  for (d in 1:Dsets) {
    log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] = log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
    log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)]   = log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
    log_mean_cup[cbegin[d]:(cbegin[d+1]-1)]    = log_mean_cup[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
    log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)]  = log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
  }
}
