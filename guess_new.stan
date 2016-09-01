////
// Cut-site normalization model: initial guess
////
functions {
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
      //print(ctr, " left=", left, " (",x[left], ") right=", right, " (", x[right], ") pos=", pos, " (", x[pos], ")");
      if ( (x[right]-z)*(x[pos]-z) > 0 ) {right = pos;} else {left = pos;}
      if (right-left==1 && x[right] > z && x[left] <= z) return left;
      pos = (right+left)/2;
    }
    return pos;
  }
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
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(SD)] XDrow_w;
  int XDrow_v[nnz(SD)];
  int XDrow_u[SD+1];
  row_vector[Krow] prowD[Dsets];
  //rescale lambdas
  real lnu;
  real ldelta;

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
  // design matrix and projector for sparse cubic spline
  //needs declaration of the following
  //int S
  //int Krow
  //vector[S] cutsites
  //int splinedegree()
  //vector[nnz(S)] Xrow_w;
  //int Xrow_v[nnz(S)];
  //int Xrow_u[S+1];
  //vector[S] row_weights;
  //row_vector[Krow] prow;
  //and row_weights needs to be filled

  ////bias spline, sparse (nu and delta have the same design)
  //BEGIN sparse calculation
  //cannot write function that modifies its arguments so we put it here
  //input: vector[S] cutsites, int Krow, int splinedegree()
  //output: vector[nnz(S)] Xrow_w, int Xrow_v[nnz(S)], int Xrow_u[S+1]
  {
    real dx; //interval length
    row_vector[Krow] t; //Krownot locations (except last)
    int x_n[Krow-splinedegree()+1]; //cumulative histogram of cutsites values
    int idx_u; //counters for filling of sparse matrix
    int idx_w;
    //
    dx = 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree()); //make it slightly larger
    t = min(cutsites) - dx*0.01 + dx * range(-splinedegree(),Krow-splinedegree()-1)';
    //get indices of cutsites which cross to the next segment and build x_n.
    x_n[1] = 1;
    x_n[2:(Krow-splinedegree())] = cumulative_hist(cutsites, t[(splinedegree()+2):]);
    x_n[Krow-splinedegree()+1] = rows(cutsites)+1;
    //
    //build spline, interval per interval
    idx_u = 1;
    idx_w = 1;
    for (i in 1:(Krow-splinedegree())) {
      //at any cutsites, there are splinedegree()+1 nonzero b-splines. Compute them.
      int xbegin;
      int xend;
      xbegin = x_n[i];
      xend = x_n[i+1]-1;
      {
        matrix[xend-xbegin+1,splinedegree()+1] tmp;
        tmp = bspl_gen(cutsites[xbegin:xend], dx, t[i:(i+splinedegree())], splinedegree());
        for (ti in 1:(xend-xbegin+1)) {
          Xrow_u[idx_u] = idx_w;
          idx_u = idx_u + 1;
          for (tj in 1:(splinedegree()+1)) {
            Xrow_w[idx_w] = tmp[ti,tj];
            Xrow_v[idx_w] = i+tj-1;
            idx_w = idx_w + 1;
          }
        }
      }
    }
    Xrow_u[idx_u] = idx_w; 
  }
  //END sparse calculation

  //projector for biases (GAM)
  row_weights = row_weights/mean(row_weights);
  prow = vector_times_csr_matrix(Krow, row_weights, Xrow_w, Xrow_v, Xrow_u);
  prow = prow / (prow * rep_vector(1,Krow));
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
    lnu = lambda_nu*tmp;
    ldelta = lambda_delta*tmp;
  }
}
parameters {
  //exposures
  real eC[Dsets];
  real eRJ[Dsets];
  real eDE[Dsets];
  //spline parameters
  vector[Krow-1] beta_nu[Biases];
  vector[Krow-1] beta_delta[Biases];
  //dispersion
  real<lower=0> alpha[Dsets];
}
transformed parameters {
  //nu
  vector[Krow] beta_nu_aug[Dsets];
  vector[Dsets*Krow] beta_nu_centered;
  vector[SD] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff[Dsets]; //2nd order difference on beta_nu_aug
  //delta
  vector[SD] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff[Dsets]; //2nd order difference on beta_delta_aug
  //means
  vector[SD] log_mean_DL;
  vector[SD] log_mean_DR;
  vector[SD] log_mean_RJ;
  //
  vector[SD] log_mean_cleft;
  vector[SD] log_mean_cright;

  //nu
  {
    for (d in 1:Dsets) {
      vector[Krow] tmp;
      beta_nu_aug[d,1] = 0;
      beta_nu_aug[d,2:] = beta_nu[XB[d]];
      tmp = beta_nu_aug[d] - (prowD[d] * beta_nu_aug[d]) * rep_vector(1,Krow);
      beta_nu_centered[((d-1)*Krow+1):(d*Krow)] = tmp;
      beta_nu_diff[d] = tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_nu = csr_matrix_times_vector(SD, Dsets*Krow, XDrow_w, XDrow_v, XDrow_u, beta_nu_centered);
  }
  
  //delta
  {
    vector[Dsets*Krow] beta_delta_centered;
    for (d in 1:Dsets) {
      vector[Krow] beta_delta_aug;
      vector[Krow] tmp;
      beta_delta_aug[1] = 0;
      beta_delta_aug[2:] = beta_delta[XB[d]];
      tmp = beta_delta_aug - (prowD[d] * beta_delta_aug) * rep_vector(1,Krow);
      beta_delta_centered[((d-1)*Krow+1):(d*Krow)] = tmp;
      beta_delta_diff[d] = tmp[:(Krow-2)]-2*tmp[2:(Krow-1)]+tmp[3:];
    }
    log_delta = csr_matrix_times_vector(SD, Dsets*Krow, XDrow_w, XDrow_v, XDrow_u, beta_delta_centered);
  }
  
  //means
  {
    //biases
    log_mean_RJ = log_nu;
    log_mean_DL = log_nu + log_delta;
    log_mean_DR = log_nu - log_delta;
    //exact counts  
    log_mean_cleft  = log_nu + log_delta;
    log_mean_cright = log_nu - log_delta;
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
  for (d in 1:Dsets) {
    rejoined[bbegin[d]:(bbegin[d+1]-1)]  ~ neg_binomial_2_log(log_mean_RJ[bbegin[d]:(bbegin[d+1]-1)], alpha[d]);
    danglingL[bbegin[d]:(bbegin[d+1]-1)] ~ neg_binomial_2_log(log_mean_DL[bbegin[d]:(bbegin[d+1]-1)], alpha[d]);
    danglingR[bbegin[d]:(bbegin[d+1]-1)] ~ neg_binomial_2_log(log_mean_DR[bbegin[d]:(bbegin[d+1]-1)], alpha[d]);
  }
  //counts
  for (d in 1:Dsets) {
    target += (bbegin[d+1]-bbegin[d]-1) * neg_binomial_2_log_lpmf(
      counts_sum_left[bbegin[d]:(bbegin[d+1]-1)] | log_mean_cleft[bbegin[d]:(bbegin[d+1]-1)], alpha[d]);
    target += (bbegin[d+1]-bbegin[d]-1) * neg_binomial_2_log_lpmf(
      counts_sum_right[bbegin[d]:(bbegin[d+1]-1)] | log_mean_cright[bbegin[d]:(bbegin[d+1]-1)], alpha[d]);
  }
  
  //// prior
  log_nu ~ cauchy(0, 1); //give high probability to [0.5:2]
  log_delta ~ cauchy(0,1);
  for (d in 1:Dsets) {
    beta_nu_diff[d] ~ normal(0,1/lnu);
    beta_delta_diff[d] ~ normal(0,1/ldelta);
    beta_nu_diff[d] ~ double_exponential(0,10/lnu);
    beta_delta_diff[d] ~ double_exponential(0,10/ldelta);
  }
} 
generated quantities {
  row_vector[Krow] genprowD[Dsets];
  genprowD <- prowD;
}