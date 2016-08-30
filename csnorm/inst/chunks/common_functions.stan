 ///////////
  //BEGIN internal use
  //
  //generate [imin, imin+1, ... imax]
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
  }
  //Compute cumulative histogram of x at left cut points q
  //i.e. index of first x value that falls in bins defined by cutpoints defined in q
  //assumes x is ordered.
  //x[hist[i]] is the first value in x that is >= q[i]
  int[] cumulative_hist(vector x, row_vector q) {
    int indices[cols(q)];
    int ix;
    int N;
    ix = 1;
    N = rows(x);
    for (iq in 1:cols(q)) {
      while (ix < N && x[ix] < q[iq]) ix = ix + 1;
      indices[iq] = ix;
    }
    return indices;
  }
  matrix bspl_gen(vector x, real dx, row_vector t, int q) {
    int N;
    int K;
    K = cols(t);
    N = rows(x);
    {
      int r[K];
      matrix[N, K] T;
      matrix[N, K] X;
      matrix[N, K] P;
      matrix[N, K] B;
      for (i in 2:K) r[i-1] = i;
      r[K] = 1;
      T = rep_matrix(t, N);
      X = rep_matrix(x, K);
      P = (X - T) / dx;
      for (i in 1:N)
        for (j in 1:K)
          B[i,j] = (T[i,j] <= X[i,j]) && (X[i,j] < T[i,j]+dx);
      for (k in 1:q)
        B = ( P .* B + (k+1-P) .* B[,r]) / k;
      return B;
    }
  }
  //END internal use
  ///////////
  
  int splinedegree() {return 3;} //set 3 for cubic spline
  
  int nnz(int N) {return N*(splinedegree()+1);} //nonzero count for design matrix
  
  //implement w^T * X  i.e. left multiply
  row_vector vector_times_csr_matrix(int K, vector w, vector Xw, int[] Xv, int[] Xu) {
    int N;
    row_vector[K] sums;
    N = rows(w);
    if (size(Xu) != N+1) reject("Xu is not of size ", N+1);
    sums = rep_row_vector(0,K);
    for (i in 1:N) {
      for (j in Xu[i]:(Xu[i+1]-1)) {
        sums[Xv[j]] = sums[Xv[j]] + Xw[j]*w[i];
      }
    }
    return sums;
  }
  
  matrix bspline(vector x, int K, int q, real xmin, real xmax) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx = 1.01*(xmax-xmin)/(K-q); //make it slightly larger
    t = xmin - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }

  int[] change_points(int[] Nvals, int levels) {
    int Nidx[levels+1];
    int i;
    int nN;
    Nidx[1] = 1;
    i = 2;
    nN = size(Nvals);
    for (N in 2:(levels+1)) {
      while (i <= nN && Nvals[i] == Nvals[i-1]) i = i + 1;
      Nidx[N] = i;
      i = i + 1;
    }
    if (Nidx[levels+1] != nN+1) reject("Nidx badly built or incorrect number of levels");
    return Nidx;
  }
  
  real neg_binomial_2_log_deviance(int[] y, vector log_mu, real alpha, vector weights) {
    vector[size(y)] y_vec;
    vector[size(y)] y_mod;
    vector[size(y)] mu;
    vector[size(y)] wt;
    //
    y_vec = to_vector(y);
    //
    if (rows(log_mu) != size(y)) {
      if (rows(log_mu) != 1) reject("size of log_mu (", rows(log_mu),
                                     ") must be either 1 or the size of y (",size(y),")");
      mu = rep_vector(exp(log_mu[1]), size(y));
    } else {
      mu = exp(log_mu);
    }
    //
    if (rows(weights) != size(y)) {
      if (rows(weights) != 1) reject("size of weights (", rows(weights),
                                     ") must be either 1 or the size of y (",size(y),")");
      wt = rep_vector(weights[1], size(y));
    } else {
      wt = weights;
    }
    //
    for (i in 1:size(y)) if (y[i]>0) {y_mod[i] = y[i];} else {y_mod[i] =1;}
    return 2*sum( ( (y_vec+alpha) .* log( (mu+alpha) ./ (y_vec+alpha) )
                    + y_vec .* log(y_mod ./ mu) ) .* wt );
  }

  // column sumn of a csr matrix
  row_vector column_sums(int K, vector Xw, int[] Xv) {
    row_vector[K] sums;
    sums = rep_row_vector(0,K);
    for (i in 1:rows(Xw)) {
      sums[Xv[i]] = sums[Xv[i]] + Xw[i];
    }
    return sums;
  }
  
  //find the index in x such as x[i]<=z and x[i+1]>z
  //startpos is a suggested starting position for i
  int bisect(real z, int startpos, vector x) {
    int pos;
    int left;
    int right;
    int nmax;
    int ctr;
    left = 1;
    right = rows(x);
    ctr = 1;
    nmax = rows(x);
    pos = startpos;
    //return edges if queried
    if (z<x[1]) return 1;
    if (z>x[nmax]) return nmax;
    //shortcut if z has same pos or the one after
    if (pos < right && x[pos] <= z && x[pos+1] > z) return pos;
    if (pos < right-1 && x[pos+1] <= z && x[pos+2] > z) return pos;
    //bisection algorithm
    while (ctr <= nmax && left != right) {
      ctr = ctr + 1;
      print(ctr, " left=", left, " (",x[left], ") right=", right, " (", x[right], ") pos=", pos, " (", x[pos], ")");
      if ( (x[right]-z)*(x[pos]-z) > 0 ) {right = pos;} else {left = pos;}
      if (right-left==1 && x[right] > z && x[left] <= z) return left;
      pos = (right+left)/2;
    }
    return pos;
  }
