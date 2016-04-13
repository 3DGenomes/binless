////
// Cut-site normalization model: fit available data
////
functions {
  ///////////
  //BEGIN internal use
  //
  //generate [imin, imin+1, ... imax]
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
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
  //parameters
  real eC;  //exposure for counts
  real eRJ; //exposure for rejoined ends
  real eDE; //exposure for dangling ends
  vector[S] log_nu; //nu_i
  vector[S] log_delta; //delta_i
  real<lower=0> lambda_diag;
}
transformed data {
  //diagonal SCAM spline, dense
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag-1] pdiag;
  
  //diagonal SCAM spline, dense
  {
    //can't do abs() on a vector so be inventive
    vector[N] tmp;
    tmp <- cutsites[cidx[2]]-cutsites[cidx[1]];
    Xdiag <- bspline(0.5*log(tmp .* tmp), Kdiag, splinedegree());
  }
  
  //projector for diagonal (SCAM)
  {
    row_vector[Kdiag] tmp;
    tmp <- rep_row_vector(1,N) * Xdiag;
    pdiag <- -tmp[2:] / (tmp * rep_vector(1,Kdiag));
  }
}
parameters {
  positive_ordered[Kdiag-1] beta_diag;
  real<lower=0> alpha;
}
transformed parameters {
  //diag
  vector[N] log_decay;
  vector[Kdiag-2] beta_diag_diff;
  //means
  vector[N] log_deltai;
  vector[N] log_deltaj;
  vector[N] base_count;
  //
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cclose;
  
  //decay
  {
    vector[Kdiag] beta_diag_centered;
    real epsilon;
    real val;
    epsilon <- -1; //+1 for increasing, -1 for decreasing spline
    val <- epsilon*pdiag*beta_diag;
    beta_diag_centered[1] <- val;
    beta_diag_centered[2:] <- val+epsilon*beta_diag;
    log_decay <- Xdiag * beta_diag_centered;
    beta_diag_diff <- beta_diag_centered[:(Kdiag-2)]-2*beta_diag_centered[2:(Kdiag-1)]+beta_diag_centered[3:];
  }

  //means
  {
    base_count <- eC + log_decay + log_nu[cidx[1]]/2 + log_nu[cidx[2]]/2;
    log_deltai <- log_delta[cidx[1]]/2;
    log_deltaj <- log_delta[cidx[2]]/2;
    log_mean_cclose <- base_count - log_deltai + log_deltaj;
    log_mean_cfar   <- base_count + log_deltai - log_deltaj;
    log_mean_cup    <- base_count + log_deltai + log_deltaj;
    log_mean_cdown  <- base_count - log_deltai - log_deltaj;
  }
}
model {
  //// likelihoods
  //counts: Close, Far, Up, Down
  counts[1] ~ neg_binomial_2_log(log_mean_cclose, alpha); // Close
  counts[2] ~ neg_binomial_2_log(log_mean_cfar, alpha); // Far
  counts[3] ~ neg_binomial_2_log(log_mean_cup, alpha); // Up
  counts[4] ~ neg_binomial_2_log(log_mean_cdown, alpha); // Down
  
  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_diag_diff ~ normal(0, 1./(alpha*lambda_diag));
}
generated quantities {
  real deviance;
  real deviance_null;
  real deviance_proportion_explained;
  deviance <- neg_binomial_2_log_deviance(counts[1], log_mean_cclose, alpha) +
              neg_binomial_2_log_deviance(counts[2], log_mean_cfar, alpha) +
              neg_binomial_2_log_deviance(counts[3], log_mean_cup, alpha) +
              neg_binomial_2_log_deviance(counts[4], log_mean_cdown, alpha);
  {
    vector[N] offset;
    offset <- rep_vector(eC, N);
    deviance_null <- neg_binomial_2_log_deviance(counts[1], offset, alpha) +
                     neg_binomial_2_log_deviance(counts[2], offset, alpha) +
                     neg_binomial_2_log_deviance(counts[3], offset, alpha) +
                     neg_binomial_2_log_deviance(counts[4], offset, alpha);
  }
  deviance_proportion_explained <- 100*(deviance_null - deviance)/deviance_null;
}
