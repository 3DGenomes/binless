////
// Cut-site normalization model: generate iota and rho given their spline parameters
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //genomic biases
  int<lower=4> Krow; //number of functions in spline base for row biases
  int<lower=1> S; //number of samples along genome
  int<lower=1> begin; //position of first and last cut site
  int<lower=begin+1> end;
  //fitted things
  vector[Krow-1] beta_iota;
  vector[Krow-1] beta_rho;
}
transformed data {
  vector[S] cutsites; //x positions, has to be within bounds of xrange
  //bias spline, sparse (iota and rho have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;
  
  ////bias spline, sparse
  cutsites = begin + range(0,S-1)*(end-begin)/(S-1);
  row_weights = rep_vector(1, S);
  #compute design matrix and projector
  #include "sparse_spline_construction.stan"
}
parameters {}
model {}
generated quantities {
  vector[S] pos;
  vector[S] log_iota;
  vector[S] log_rho;
  
  pos = cutsites;
  //iota
  {
    vector[Krow] beta_iota_aug;
    vector[Krow] beta_iota_centered;
    beta_iota_aug[1] = 0;
    beta_iota_aug[2:] = beta_iota;
    beta_iota_centered = beta_iota_aug - (prow * beta_iota_aug) * rep_vector(1,Krow);
    log_iota = csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_iota_centered);
  }
  //rho
  {
    vector[Krow] beta_rho_aug;
    vector[Krow] beta_rho_centered;
    beta_rho_aug[1] = 0;
    beta_rho_aug[2:] = beta_rho;
    beta_rho_centered = beta_rho_aug - (prow * beta_rho_aug) * rep_vector(1,Krow);
    log_rho = csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_rho_centered);
  }
}
