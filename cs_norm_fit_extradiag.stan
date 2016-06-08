////
// Cut-site normalization model: only fit diagonal decay and eC, given nu and delta
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
  int<lower=1> S1; //number of cut sites on x and y axes
  int<lower=1> S2;
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
  //fitted parameters
  vector[S1] log_nu1; // log(nu)
  vector[S2] log_nu2;
  vector[S1] log_delta1; // log(delta)
  vector[S2] log_delta2;
}
transformed data {
  //diagonal SCAM spline, dense, exact
  matrix[Nclose+Nfar+Nup+Ndown,Kdiag] Xdiag;
  row_vector[Kdiag] pdiag;
  vector[Nclose+Nfar+Nup+Ndown] diag_weights;
  
  //diagonal SCAM spline, dense, exact
  {
    vector[Nclose+Nfar+Nup+Ndown] tmp;
    tmp[:Nclose] <- dist_close;
    tmp[(Nclose+1):(Nclose+Nfar)] <- dist_far;
    tmp[(Nclose+Nfar+1):(Nclose+Nfar+Nup)] <- dist_up;
    tmp[(Nclose+Nfar+Nup+1):] <- dist_down;
    Xdiag <- bspline(log(tmp), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    diag_weights <- rep_vector(1, Nclose+Nfar+Nup+Ndown);
    diag_weights <- diag_weights/mean(diag_weights);
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
  vector[Nclose] log_decay_close;
  vector[Nfar] log_decay_far;
  vector[Nup] log_decay_up;
  vector[Ndown] log_decay_down;
  vector[Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff;
  //
  vector[Nclose] log_mean_cclose;
  vector[Nfar] log_mean_cfar;
  vector[Nup] log_mean_cup;
  vector[Ndown] log_mean_cdown;
  
  //decay
  {
    vector[Kdiag] beta_diag_aug;
    vector[Nclose+Nfar+Nup+Ndown] log_decay;
    real epsilon;
    epsilon <- -1; //decreasing spline
    beta_diag_aug[1] <- 0;
    beta_diag_aug[2:] <- beta_diag;
    beta_diag_centered <- epsilon * (beta_diag_aug - (pdiag * beta_diag_aug) * rep_vector(1, Kdiag));
    log_decay <- Xdiag * beta_diag_centered;
    log_decay_close <- log_decay[:Nclose];
    log_decay_far <- log_decay[(Nclose+1):(Nclose+Nfar)];
    log_decay_up <- log_decay[(Nclose+Nfar+1):(Nclose+Nfar+Nup)];
    log_decay_down <- log_decay[(Nclose+Nfar+Nup+1):];
    beta_diag_diff <- beta_diag_centered[:(Kdiag-2)]-2*beta_diag_centered[2:(Kdiag-1)]+beta_diag_centered[3:];
  }

  //means
  {
    //exact counts  
    log_mean_cclose <- eC + log_decay_close + (log_nu1 - log_delta1)[index_close[1]] + (log_nu2 + log_delta2)[index_close[2]];
    log_mean_cfar   <- eC + log_decay_far   + (log_nu1 + log_delta1)[index_far[1]]   + (log_nu2 - log_delta2)[index_far[2]];
    log_mean_cup    <- eC + log_decay_up    + (log_nu1 + log_delta1)[index_up[1]]    + (log_nu2 + log_delta2)[index_up[2]];
    log_mean_cdown  <- eC + log_decay_down  + (log_nu1 - log_delta1)[index_down[1]]  + (log_nu2 - log_delta2)[index_down[2]];
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
