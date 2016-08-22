////
// Cut-site normalization model: fit available data, partial mean-field approximation
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
  int<lower=1,upper=SD> bbegin[Dsets+1]; //bbegin[i]=j: dataset i starts at j
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
  int<lower=1,upper=N> cbegin[Dsets+1]; //cbegin[i]=j: dataset i starts at j
  int<lower=1,upper=S> cidx[2,N]; //indices of rsite pairs, referring to vector[SD] cutsitesD
  vector<lower=dmin,upper=dmax>[N] dist; //genomic distance between rsites
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
  real<lower=0> weight[Dsets]; //in case of subsampling
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(SD)] XDrow_w;
  int XDrow_v[nnz(SD)];
  int XDrow_u[SD+1];
  row_vector[Krow] prowD[Dsets];
  //diagonal SCAM spline, dense
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag] pdiagD[Dsets];
  //scaling factor for genomic lambdas
  real lfac;
  
  //Bias GAM spline, sparse
  {
    vector[SD] row_weightsD;
    int nnzs[Dsets+1];
    row_weightsD <- rep_vector(3, SD); //dangling L/R + rejoined
    for (i in 1:N) {
      row_weightsD[cidx[1,i]] <- row_weightsD[cidx[1,i]] + 1;
      row_weightsD[cidx[2,i]] <- row_weightsD[cidx[2,i]] + 1;
    }
    nnzs[1] <- 1;
    for (i in 1:Dsets) nnzs[i+1] <- nnzs[i]+nnz(bbegin[i+1]-bbegin[i]);
    
    #compute design matrix and projector
    for (d in 1:Dsets) {
      int S;
      S <- bbegin[d+1]-bbegin[d];
      {
        vector[S] cutsites;
        vector[nnz(S)] Xrow_w;
        int Xrow_v[nnz(S)];
        int Xrow_u[S+1];
        vector[S] row_weights;
        row_vector[Krow] prow;
        cutsites <- cutsitesD[bbegin[d]:(bbegin[d+1]-1)];
        row_weights <- row_weightsD[bbegin[d]:(bbegin[d+1]-1)];
        #include "sparse_spline_construction.stan"
        XDrow_w[nnzs[d]:(nnzs[d+1]-1)] <- Xrow_w;
        XDrow_v[nnzs[d]:(nnzs[d+1]-1)] <- Xrow_v + Krow*(d-1);
        XDrow_u[bbegin[d]:bbegin[d+1]] <- Xrow_u + bbegin[d]-1;
        prowD[d] <- prow;
      }
    }
  }
  //diagonal SCAM spline, dense
  {
    Xdiag <- bspline(log(dist), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    for (d in 1:Dsets) {
      int sz;
      sz <- cbegin[d+1]-cbegin[d];
      {
        row_vector[sz] diag_weights;
        vector[Kdiag] pdiag;
        diag_weights <- rep_row_vector(1,sz);
        diag_weights <- diag_weights/mean(diag_weights);
        pdiag <- diag_weights * Xdiag[cbegin[d]:(cbegin[d+1]-1),:];
        pdiagD[d] <- pdiag / (pdiag * rep_vector(1,Kdiag));
      }
    }
  }
  
  //scaling factor for genomic lambdas
  lfac <- 30000*Krow/(max(cutsites)-min(cutsites));
}
parameters {
  //exposures
  real eC[Dsets];
  real eRJ[Dsets];
  real eDE[Dsets];
  //spline parameters
  vector[Krow-1] beta_nu[Biases];
  vector[Krow-1] beta_delta[Biases];
  positive_ordered[Kdiag-1] beta_diag[Decays];
  //dispersion
  real<lower=0> alpha;
  //length scales
  real<lower=0> lambda_nu[Biases];
  real<lower=0> lambda_delta[Biases];
  real<lower=0> lambda_diag[Decays];
}
transformed parameters {
  //nu
  vector[SD] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff[Dsets]; //2nd order difference on beta_nu_aug
  //delta
  vector[SD] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff[Dsets]; //2nd order difference on beta_delta_aug
  //diag
  vector[N] log_decay;
  vector[Dsets*Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff[Dsets];
  //means
  vector[SD] log_mean_DL;
  vector[SD] log_mean_DR;
  vector[SD] log_mean_RJ;
  //
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  
  //nu
  {
    vector[Dsets*Krow] beta_nu_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_nu_aug;
      vector[Krow] tmp;
      beta_nu_aug[1] <- 0;
      beta_nu_aug[2:] <- beta_nu[XB[d]];
      tmp <- beta_nu_aug - (prowD[d] * beta_nu_aug) * rep_vector(1,Krow);
      beta_nu_centered[((d-1)*Krow+1):(d*Krow)] <- tmp;
      beta_nu_diff[d] <- tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_nu <- csr_matrix_times_vector(SD, Dsets*Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
  }
  
  //delta
  {
    vector[Dsets*Krow] beta_delta_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_delta_aug;
      vector[Krow] tmp;
      beta_delta_aug[1] <- 0;
      beta_delta_aug[2:] <- beta_delta[XB[d]];
      tmp <- beta_delta_aug - (prowD[d] * beta_delta_aug) * rep_vector(1,Krow);
      beta_delta_centered[((d-1)*Krow+1):(d*Krow)] <- tmp;
      beta_delta_diff[d] <- tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_delta <- csr_matrix_times_vector(SD, Dsets*Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
  }
  
  //decay
  for (d in 1:Dsets) {
    vector[Kdiag] beta_diag_aug;
    vector[Kdiag] tmp;
    real epsilon;
    epsilon <- -1; //decreasing spline
    beta_diag_aug[1] <- 0;
    beta_diag_aug[2:] <- beta_diag[XD[d]];
    tmp <- epsilon * (beta_diag_aug - (pdiagD[d] * beta_diag_aug) * rep_vector(1, Kdiag));
    beta_diag_centered[((d-1)*Kdiag+1):(d*Kdiag)] <- tmp;
    beta_diag_diff[d] <- tmp[:(Kdiag-2)]-2*tmp[2:(Kdiag-1)]+tmp[3:];
  }
  log_decay <- Xdiag * beta_diag_centered;
    
  //means
  {
    //biases
    log_mean_RJ <- log_nu;
    log_mean_DL <- log_nu + log_delta;
    log_mean_DR <- log_nu - log_delta;
    //exact counts  
    log_mean_cclose <- log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cfar   <- log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
    log_mean_cup    <- log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cdown  <- log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
    //add exposures
    for (d in 1:Dsets) {
      log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)]     <- log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)] + eRJ[d];
      log_mean_DL[bbegin[d]:(bbegin[d+1]-1)]     <- log_mean_DL[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_DR[bbegin[d]:(bbegin[d+1]-1)]     <- log_mean_DR[bbegin[d]:(bbegin[d+1]-1)] + eDE[d];
      log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] <- log_mean_cclose[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
      log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)]   <- log_mean_cfar[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
      log_mean_cup[cbegin[d]:(cbegin[d+1]-1)]    <- log_mean_cup[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
      log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)]  <- log_mean_cdown[cbegin[d]:(cbegin[d+1]-1)] + eC[d];
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
  increment_log_prob(weight*neg_binomial_2_log_log(counts_close, log_mean_cclose, alpha));
  increment_log_prob(weight*neg_binomial_2_log_log(counts_far, log_mean_cfar, alpha));
  increment_log_prob(weight*neg_binomial_2_log_log(counts_up, log_mean_cup, alpha));
  increment_log_prob(weight*neg_binomial_2_log_log(counts_down, log_mean_cdown, alpha));
  
  //// Priors
  for (d in 1:Dsets) {
    //P-spline prior on the differences (K-2 params)
    //warning on jacobian can be ignored
    //see GAM, Wood (2006), section 4.8.2 (p.187)
    beta_nu_diff[d] ~ normal(0, 1/(lfac*lambda_nu[XB[d]]));
    beta_delta_diff[d] ~ normal(0, 1/(lfac*lambda_delta[XB[d]]));
    beta_diag_diff[d] ~ normal(0, 1/lambda_diag[XD[d]]);
    //weak lasso on bias splines (non-bayesian)
    beta_nu_diff[d] ~ double_exponential(0,10/(lfac*lambda_nu[XB[d]]));
    beta_delta_diff[d] ~ double_exponential(0,10/(lfac*lambda_delta[XB[d]]));
  }
  //cauchy hyperprior
  lambda_nu ~ cauchy(0,0.1);
  lambda_delta ~ cauchy(0,0.1);
  lambda_diag ~ cauchy(0,0.1);
}
