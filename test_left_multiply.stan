functions {
  //implement w^T * X  i.e. left multiply
  row_vector vector_times_csr_matrix(int K, vector w, vector Xw, int[] Xv, int[] Xu) {
    int N;
    row_vector[K] sums;
    N <- rows(w);
    if (size(Xu) != N+1) reject("Xu is not of size ", N+1);
    sums <- rep_row_vector(0,K);
    for (i in 1:N) {
      for (j in Xu[i]:(Xu[i+1]-1)) {
        sums[Xv[j]] <- sums[Xv[j]] + Xw[j]*w[i];
      }
    }
    return sums;
  }
}
////////////////////////////////////////////////////////////////
data {
  int K;
  int N;
  int nnz;
  matrix[N,K] densemat;
  vector[N] weights;
}
parameters {}
model {}
generated quantities {
  row_vector[K] densevec;
  row_vector[K] sparsevec;
  real diff;
  densevec <- weights' * densemat;
  sparsevec <- vector_times_csr_matrix(K, weights, csr_extract_w(densemat), csr_extract_v(densemat),
                                       csr_extract_u(densemat));
  diff <- sum(densevec-sparsevec);
}
