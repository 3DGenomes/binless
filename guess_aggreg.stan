////
// Cut-site normalization model: fit available data, partial mean-field approximation
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
    ix <- 1;
    N <- rows(x);
    for (iq in 1:cols(q)) {
      while (ix < N && x[ix] < q[iq]) ix <- ix + 1;
      indices[iq] <- ix;
    }
    return indices;
  }
  matrix bspl_gen(vector x, real dx, row_vector t, int q) {
    int N;
    int K;
    K <- cols(t);
    N <- rows(x);
    {
      int r[K];
      matrix[N, K] T;
      matrix[N, K] X;
      matrix[N, K] P;
      matrix[N, K] B;
      for (i in 2:K) r[i-1] <- i;
      r[K] <- 1;
      T <- rep_matrix(t, N);
      X <- rep_matrix(x, K);
      P <- (X - T) / dx;
      for (i in 1:N)
        for (j in 1:K)
          B[i,j] <- (T[i,j] <= X[i,j]) && (X[i,j] < T[i,j]+dx);
      for (k in 1:q)
        B <- ( P .* B + (k+1-P) .* B[,r]) / k;
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
  
  matrix bspline(vector x, int K, int q, real xmin, real xmax) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(xmax-xmin)/(K-q); //make it slightly larger
    t <- xmin - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }

  int[] change_points(int[] Nvals, int levels) {
    int Nidx[levels+1];
    int i;
    int nN;
    Nidx[1] <- 1;
    i <- 2;
    nN <- size(Nvals);
    for (N in 2:(levels+1)) {
      while (i <= nN && Nvals[i] == Nvals[i-1]) i <- i + 1;
      Nidx[N] <- i;
      i <- i + 1;
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
    y_vec <- to_vector(y);
    //
    if (rows(log_mu) != size(y)) {
      if (rows(log_mu) != 1) reject("size of log_mu (", rows(log_mu),
                                     ") must be either 1 or the size of y (",size(y),")");
      mu <- rep_vector(exp(log_mu[1]), size(y));
    } else {
      mu <- exp(log_mu);
    }
    //
    if (rows(weights) != size(y)) {
      if (rows(weights) != 1) reject("size of weights (", rows(weights),
                                     ") must be either 1 or the size of y (",size(y),")");
      wt <- rep_vector(weights[1], size(y));
    } else {
      wt <- weights;
    }
    //
    for (i in 1:size(y)) if (y[i]>0) {y_mod[i] <- y[i];} else {y_mod[i] <-1;}
    return 2*sum( ( (y_vec+alpha) .* log( (mu+alpha) ./ (y_vec+alpha) )
                    + y_vec .* log(y_mod ./ mu) ) .* wt );
  }

  // column sumn of a csr matrix
  row_vector column_sums(int K, vector Xw, int[] Xv) {
    row_vector[K] sums;
    sums <- rep_row_vector(0,K);
    for (i in 1:rows(Xw)) {
      sums[Xv[i]] <- sums[Xv[i]] + Xw[i];
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
    left <- 1;
    right <- rows(x);
    ctr <- 1;
    nmax <- rows(x);
    pos <- startpos;
    //return edges if queried
    if (z<x[1]) return 1;
    if (z>x[nmax]) return nmax;
    //shortcut if z has same pos or the one after
    if (pos < right && x[pos] <= z && x[pos+1] > z) return pos;
    if (pos < right-1 && x[pos+1] <= z && x[pos+2] > z) return pos;
    //bisection algorithm
    while (ctr <= nmax && left != right) {
      ctr <- ctr + 1;
      print(ctr, " left=", left, " (",x[left], ") right=", right, " (", x[right], ") pos=", pos, " (", x[pos], ")");
      if ( (x[right]-z)*(x[pos]-z) > 0 ) {right <- pos;} else {left <- pos;}
      if (right-left==1 && x[right] > z && x[left] <= z) return left;
      pos <- (right+left)/2;
    }
    return pos;
  }}
////////////////////////////////////////////////////////////////
data {
  //genomic biases
  int<lower=4> Krow; //number of functions in spline base for row biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  int<lower=0> rejoined[S];
  int<lower=0> danglingL[S];
  int<lower=0> danglingR[S];
  //count sums
  int<lower=0> counts_sum_left[S];
  int<lower=0> counts_sum_right[S];
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;

  row_weights <- rep_vector(3, S); //dangling L/R + rejoined
  #compute design matrix and projector
// design matrix and projector for sparse cubic spline
  //needs declaration of the following
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
    dx <- 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree()); //make it slightly larger
    t <- min(cutsites) - dx*0.01 + dx * range(-splinedegree(),Krow-splinedegree()-1)';
    //get indices of cutsites which cross to the next segment and build x_n.
    x_n[1] <- 1;
    x_n[2:(Krow-splinedegree())] <- cumulative_hist(cutsites, t[(splinedegree()+2):]);
    x_n[Krow-splinedegree()+1] <- rows(cutsites)+1;
    //
    //build spline, interval per interval
    idx_u <- 1;
    idx_w <- 1;
    for (i in 1:(Krow-splinedegree())) {
      //at any cutsites, there are splinedegree()+1 nonzero b-splines. Compute them.
      int xbegin;
      int xend;
      xbegin <- x_n[i];
      xend <- x_n[i+1]-1;
      {
        matrix[xend-xbegin+1,splinedegree()+1] tmp;
        tmp <- bspl_gen(cutsites[xbegin:xend], dx, t[i:(i+splinedegree())], splinedegree());
        for (ti in 1:(xend-xbegin+1)) {
          Xrow_u[idx_u] <- idx_w;
          idx_u <- idx_u + 1;
          for (tj in 1:(splinedegree()+1)) {
            Xrow_w[idx_w] <- tmp[ti,tj];
            Xrow_v[idx_w] <- i+tj-1;
            idx_w <- idx_w + 1;
          }
        }
      }
    }
    Xrow_u[idx_u] <- idx_w; 
  }
  //END sparse calculation

  //projector for biases (GAM)
  row_weights <- row_weights/mean(row_weights);
  prow <- vector_times_csr_matrix(Krow, row_weights, Xrow_w, Xrow_v, Xrow_u);
  prow <- prow / (prow * rep_vector(1,Krow));}
parameters {
  //exposures
  real eC;
  real eRJ;
  real eDE;
  //spline parameters
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  //dispersion
  real<lower=0> alpha_bias;
  real<lower=0> alpha_count;
  //length scales
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
}
transformed parameters {
  //nu
  vector[S] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff; //2nd order difference on beta_nu_aug
  //delta
  vector[S] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff; //2nd order difference on beta_delta_aug
  //means
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  vector[S] log_mean_RJ;
  //
  vector[S] log_mean_cleft;
  vector[S] log_mean_cright;

  //nu
  {
    vector[Krow] beta_nu_aug;
    vector[Krow] beta_nu_centered;
    beta_nu_aug[1] <- 0;
    beta_nu_aug[2:] <- beta_nu;
    beta_nu_centered <- beta_nu_aug - (prow * beta_nu_aug) * rep_vector(1,Krow);
    log_nu <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
    beta_nu_diff <- beta_nu_centered[:(Krow-2)]-2*beta_nu_centered[2:(Krow-1)]+beta_nu_centered[3:];
  }
  //delta
  {
    vector[Krow] beta_delta_aug;
    vector[Krow] beta_delta_centered;
    beta_delta_aug[1] <- 0;
    beta_delta_aug[2:] <- beta_delta;
    beta_delta_centered <- beta_delta_aug - (prow * beta_delta_aug) * rep_vector(1,Krow);
    log_delta <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
    beta_delta_diff <- beta_delta_centered[:(Krow-2)]-2*beta_delta_centered[2:(Krow-1)]+beta_delta_centered[3:];
  }

  //means
  {
    //biases
    log_mean_RJ <- log_nu + eRJ;
    log_mean_DL <- log_nu + eDE + log_delta;
    log_mean_DR <- log_nu + eDE - log_delta;
    //exact counts  
    log_mean_cleft  <- log(2*S) + eC + log_nu + log_delta;
    log_mean_cright <- log(2*S) + eC + log_nu - log_delta;
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha_bias);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha_bias);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha_bias);
  
  //counts
  counts_sum_left  ~ neg_binomial_2_log(log_mean_cleft, alpha_count);
  counts_sum_right ~ neg_binomial_2_log(log_mean_cright, alpha_count);
  
  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_nu_diff ~ normal(0, 1./lambda_nu);
  beta_delta_diff ~ normal(0, 1./lambda_delta);
}
