////
// Cut-site normalization model: predict mean for each provided count
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
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Kdiag; //number of functions in spline base for diagonal decay
  //biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  //distance bounds
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //counts : explicit
  int<lower=0> Nclose; //number of close counts modelled explicitly
  int<lower=0> counts_close[Nclose]; //value of the count
  int<lower=0> index_close[2,Nclose]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Nclose] dist_close; //genomic distance between rsites
  //
  int<lower=0> Nfar; //number of far counts modelled explicitly
  int<lower=0> counts_far[Nfar]; //value of the count
  int<lower=0> index_far[2,Nfar]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Nfar] dist_far; //genomic distance between rsites
  //
  int<lower=0> Nup; //number of upstream counts modelled explicitly
  int<lower=0> counts_up[Nup]; //value of the count
  int<lower=0> index_up[2,Nup]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Nup] dist_up; //genomic distance between rsites
  //
  int<lower=0> Ndown; //number of downstream counts modelled explicitly
  int<lower=0> counts_down[Ndown]; //value of the count
  int<lower=0> index_down[2,Ndown]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[Ndown] dist_down; //genomic distance between rsites
  //
  //estimated parameters
  real eC;  //exposure for counts
  vector[S] log_nu; //take nu and delta directly to avoid base reconstruction
  vector[S] log_delta; 
  vector[Kdiag] beta_diag_centered; //need to build spline base
  real<lower=0> alpha;
}
parameters {}
model {}
generated quantities {
  //diag
  vector[Nclose] log_decay_close;
  vector[Nfar] log_decay_far;
  vector[Nup] log_decay_up;
  vector[Ndown] log_decay_down;
  //means
  vector[Nclose] log_mean_cclose;
  vector[Nfar] log_mean_cfar;
  vector[Nup] log_mean_cup;
  vector[Ndown] log_mean_cdown;
  
  //decay: diagonal SCAM spline, dense
  {
    matrix[Nclose+Nfar+Nup+Ndown,Kdiag] Xdiag;
    vector[Nclose+Nfar+Nup+Ndown] log_decay;
    //design matrix
    {
      vector[Nclose+Nfar+Nup+Ndown] tmpN;
      tmpN[:Nclose] <- dist_close;
      tmpN[(Nclose+1):(Nclose+Nfar)] <- dist_far;
      tmpN[(Nclose+Nfar+1):(Nclose+Nfar+Nup)] <- dist_up;
      tmpN[(Nclose+Nfar+Nup+1):] <- dist_down;
      Xdiag <- bspline(log(tmpN), Kdiag, splinedegree(), log(dmin), log(dmax));
    }
    //decay
    log_decay <- Xdiag * beta_diag_centered;
    log_decay_close <- log_decay[:Nclose];
    log_decay_far <- log_decay[(Nclose+1):(Nclose+Nfar)];
    log_decay_up <- log_decay[(Nclose+Nfar+1):(Nclose+Nfar+Nup)];
    log_decay_down <- log_decay[(Nclose+Nfar+Nup+1):];
  }
  
  //means
  log_mean_cclose <- eC + log_decay_close + (log_nu - log_delta)[index_close[1]] 
                                         + (log_nu + log_delta)[index_close[2]];
  log_mean_cfar   <- eC + log_decay_far   + (log_nu + log_delta)[index_far[1]]   
                                          + (log_nu - log_delta)[index_far[2]];
  log_mean_cup    <- eC + log_decay_up    + (log_nu + log_delta)[index_up[1]]    
                                          + (log_nu + log_delta)[index_up[2]];
  log_mean_cdown  <- eC + log_decay_down  + (log_nu - log_delta)[index_down[1]]  
                                          + (log_nu - log_delta)[index_down[2]];
}
