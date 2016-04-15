////
// Cut-site normalization model: predict for full dataset
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
  
  matrix bspline(vector x, int K, int q) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(max(x)-min(x))/(K-q); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }

  real neg_binomial_2_log_deviance(int[] y, vector log_mu, real alpha) {
    vector[size(y)] y_vec;
    vector[size(y)] y_mod;
    vector[rows(log_mu)] mu;
    if (rows(log_mu) != size(y)) reject("sizes of y (",size(y),") and log_mu (", rows(log_mu),
                                        ") must match")
    y_vec <- to_vector(y);
    for (i in 1:size(y)) if (y[i]>0) {y_mod[i] <- y[i];} else {y_mod[i] <-1;}
    mu <- exp(log_mu);
    return 2*sum( (y_vec+alpha) .* log( (mu+alpha) ./ (y_vec+alpha) )
                  + y_vec .* log(y_mod ./ mu));
  }
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Kdiag; //number of functions in spline base for diagonal decay
  //biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  //counts
  int<lower=1> N; //number of data points
  int<lower=0> counts[4,N]; //raw counts: Close, Far, Up, Down
  int<lower=1,upper=S> cidx[2,N]; //index of its associated cut site
  //estimated parameters
  real eC;  //exposure for counts
  vector[S] log_nu; //take nu and delta directly to avoid base reconstruction
  vector[S] log_delta; 
  vector[Kdiag-1] beta_diag; //need to build spline base
                             //downcast to vector because positive_ordered is too strict
  real<lower=0> alpha;
}
transformed data {
  //nu
  vector[N] log_nui;
  vector[N] log_nuj;
  //delta
  vector[N] log_deltai;
  vector[N] log_deltaj;

   //nu
  log_nui <- log_nu[cidx[1]];
  log_nuj <- log_nu[cidx[2]];
  //delta
  log_deltai <- log_delta[cidx[1]];
  log_deltaj <- log_delta[cidx[2]];
}
parameters {}
model {}
generated quantities {
  //diag
  vector[N] log_decay;
  //means
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cclose;
  
  //decay: diagonal SCAM spline, dense
  {
    matrix[N,Kdiag] Xdiag;
    row_vector[Kdiag-1] pdiag;
    vector[N] tmpN;
    row_vector[Kdiag] tmpK;
    vector[Kdiag] beta_diag_centered;
    real epsilon;
    real val;
    //X: can't do abs() on a vector so be inventive
    tmpN <- cutsites[cidx[2]]-cutsites[cidx[1]];
    Xdiag <- bspline(0.5*log(tmpN .* tmpN), Kdiag, splinedegree());
    //projector for diagonal (SCAM)
    tmpK <- rep_row_vector(1,N) * Xdiag;
    pdiag <- -tmpK[2:] / (tmpK * rep_vector(1,Kdiag));
    //decay
    epsilon <- -1; //+1 for increasing, -1 for decreasing spline
    val <- epsilon*pdiag*beta_diag;
    beta_diag_centered[1] <- val;
    beta_diag_centered[2:] <- val+epsilon*beta_diag;
    log_decay <- Xdiag * beta_diag_centered;
  }

  //means
  {
    vector[N] base_count;
    base_count <- eC + log_decay + log_nui + log_nuj;
    log_mean_cclose <- base_count - log_deltai + log_deltaj;
    log_mean_cfar   <- base_count + log_deltai - log_deltaj;
    log_mean_cup    <- base_count + log_deltai + log_deltaj;
    log_mean_cdown  <- base_count - log_deltai - log_deltaj;
  }
}
