////
// Cut-site normalization model: only fit diagonal decay and eC, given nu and delta
////
functions {
  ///////////
  //BEGIN internal use
  //
  #include "range.stan"
  //
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
  //genomic biases
  int<lower=1> S1; //number of cut sites on x and y axes
  int<lower=1> S2;
  //decay bias
  int Kdiag; //number of functions in spline base for diagonal decay
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //counts : explicit
  int<lower=0> N; //number of counts
  int<lower=1> cidx[2,N]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[N] dist; //genomic distance between rsites
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
  //fitted parameters
  vector[S1] log_nu1; // log(nu)
  vector[S2] log_nu2;
  vector[S1] log_delta1; // log(delta)
  vector[S2] log_delta2;
}
transformed data {
  //diagonal SCAM spline, dense, exact
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag] pdiag;
  vector[N] diag_weights;
  if (max(cidx[1])>S1) {reject("first index larger than rsites");}
  if (max(cidx[2])>S2) {reject("second index larger than rsites");}
  //diagonal SCAM spline, dense, exact
  {
    Xdiag <- bspline(log(dist), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    diag_weights <- rep_vector(1, N);
    pdiag <- diag_weights' * Xdiag;
    pdiag <- pdiag / (pdiag * rep_vector(1,Kdiag));
  }
}
parameters {
  //exposures
  real eC;
  positive_ordered[Kdiag-1] beta_diag;
  //deviance
  real<lower=0> alpha;
  //length scales
  real<lower=0> lambda_diag;
}
transformed parameters {
  //diag
  vector[N] log_decay;
  vector[Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff;
  //
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  
  //decay
  {
    vector[Kdiag] beta_diag_aug;
    real epsilon;
    epsilon <- -1; //decreasing spline
    beta_diag_aug[1] <- 0;
    beta_diag_aug[2:] <- beta_diag;
    beta_diag_centered <- epsilon * (beta_diag_aug - (pdiag * beta_diag_aug) * rep_vector(1, Kdiag));
    log_decay <- Xdiag * beta_diag_centered;
    beta_diag_diff <- beta_diag_centered[:(Kdiag-2)]-2*beta_diag_centered[2:(Kdiag-1)]+beta_diag_centered[3:];
  }

  //means
  {
    //exact counts  
    log_mean_cclose <- eC + log_decay + (log_nu1 - log_delta1)[cidx[1]] + (log_nu2 + log_delta2)[cidx[2]];
    log_mean_cfar   <- eC + log_decay + (log_nu1 + log_delta1)[cidx[1]] + (log_nu2 - log_delta2)[cidx[2]];
    log_mean_cup    <- eC + log_decay + (log_nu1 + log_delta1)[cidx[1]] + (log_nu2 + log_delta2)[cidx[2]];
    log_mean_cdown  <- eC + log_decay + (log_nu1 - log_delta1)[cidx[1]] + (log_nu2 - log_delta2)[cidx[2]];
  }
}
model {
  //// Exact likelihoods
  //counts: Close, Far, Up, Down
  counts_close ~ neg_binomial_2_log(log_mean_cclose, alpha); // Close
  counts_far   ~ neg_binomial_2_log(log_mean_cfar, alpha); // Far
  counts_up    ~ neg_binomial_2_log(log_mean_cup, alpha); // Up
  counts_down  ~ neg_binomial_2_log(log_mean_cdown, alpha); // Down
  

  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_diag_diff ~ normal(0, 1./lambda_diag);
}
generated quantities {
  real deviance;
  real deviance_null;
  real deviance_proportion_explained;
  #deviances
  deviance <- neg_binomial_2_log_deviance(counts_close, log_mean_cclose, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_far, log_mean_cfar, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_up, log_mean_cup, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_down, log_mean_cdown, alpha, rep_vector(1,1));
  #null deviances
  deviance_null <- neg_binomial_2_log_deviance(counts_close, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_far, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_up, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_down, rep_vector(eC,1), alpha, rep_vector(1,1));
  #proportions explained
  deviance_proportion_explained <- 100*(deviance_null - deviance)/deviance_null;
}
