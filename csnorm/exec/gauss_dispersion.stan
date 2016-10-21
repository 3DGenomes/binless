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
  int<lower=4> Krow; //number of functions in spline base for row biases
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
  vector[Krow-1] beta_iota[Biases];
  vector[Krow-1] beta_rho[Biases];
  positive_ordered[Kdiag-1] beta_diag[Decays];
}
transformed data {
  //bias spline, sparse (iota and rho have the same design)
  vector[nnz(SD)] XDrow_w;
  int XDrow_v[nnz(SD)];
  int XDrow_u[SD+1];
  row_vector[Krow] prowD[Dsets];
  //diagonal SCAM spline, dense
  matrix[N,Dsets*Kdiag] Xdiag;
  row_vector[Kdiag] pdiagD[Dsets];
  //iota
  vector[SD] log_iota; // log(iota)
  vector[Krow-2] beta_iota_diff[Dsets]; //2nd order difference on beta_iota_aug
  //rho
  vector[SD] log_rho; // log(rho)
  vector[Krow-2] beta_rho_diff[Dsets]; //2nd order difference on beta_rho_aug
  //diag
  vector[N] ldec;
  vector[Dsets*Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff[Dsets];
  
  //Bias GAM spline, sparse
  {
    vector[SD] row_weightsD;
    int nnzs[Dsets+1];
    row_weightsD = rep_vector(3, SD); //dangling L/R + rejoined
    for (i in 1:N) {
      row_weightsD[cidx[1,i]] = row_weightsD[cidx[1,i]] + 1;
      row_weightsD[cidx[2,i]] = row_weightsD[cidx[2,i]] + 1;
    }
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

