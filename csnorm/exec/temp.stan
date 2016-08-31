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
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
}
parameters {}
model {}
generated quantities {
  //bias spline, sparse (nu and delta have the same design)
  vector[4*SD] XDrow_w;
  int XDrow_v[4*SD];
  int XDrow_u[SD+1];
  //same but through dense matrix
  vector[4*SD] YDrow_w;
  int YDrow_v[4*SD];
  int YDrow_u[SD+1];
  
  //Bias GAM spline, sparse
  {
    matrix[SD, Dsets*Krow] Ydense;
    vector[SD] row_weightsD;
    int nnzs[Dsets+1];
    Ydense = rep_matrix(0, SD, Dsets*Krow);
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
        Ydense[bbegin[d]:(bbegin[d+1]-1),((d-1)*Krow+1):(d*Krow)] = csr_to_dense_matrix(S,Krow,Xrow_w,Xrow_v,Xrow_u);
        XDrow_w[nnzs[d]:(nnzs[d+1]-1)] = Xrow_w;
        for (i in 1:size(Xrow_v)) Xrow_v[i] = Xrow_v[i] + Krow*(d-1);
        XDrow_v[nnzs[d]:(nnzs[d+1]-1)] = Xrow_v;
      }
    }
    XDrow_u[1] = 1;
    for (i in 1:SD) XDrow_u[i+1] = XDrow_u[i]+nnz(1);
    YDrow_w = csr_extract_w(Ydense);
    YDrow_v = csr_extract_v(Ydense);
    YDrow_u = csr_extract_u(Ydense);
  }
}
