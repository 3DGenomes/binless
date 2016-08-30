////
// Cut-site normalization model: generate nu and delta given their spline parameters
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
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
}
transformed data {
  vector[S] cutsites; //x positions, has to be within bounds of xrange
  //bias spline, sparse (nu and delta have the same design)
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
  vector[S] log_nu;
  vector[S] log_delta;
  
  pos = cutsites;
  //nu
  {
    vector[Krow] beta_nu_aug;
    vector[Krow] beta_nu_centered;
    beta_nu_aug[1] = 0;
    beta_nu_aug[2:] = beta_nu;
    beta_nu_centered = beta_nu_aug - (prow * beta_nu_aug) * rep_vector(1,Krow);
    log_nu = csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
  }
  //delta
  {
    vector[Krow] beta_delta_aug;
    vector[Krow] beta_delta_centered;
    beta_delta_aug[1] = 0;
    beta_delta_aug[2:] = beta_delta;
    beta_delta_centered = beta_delta_aug - (prow * beta_delta_aug) * rep_vector(1,Krow);
    log_delta = csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
      }
}
