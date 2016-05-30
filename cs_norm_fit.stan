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
  
  row_vector column_sums(int K, vector Xw, int[] Xv) {
    row_vector[K] sums;
    sums <- rep_row_vector(0,K);
    for (i in 1:rows(Xw)) {
      sums[Xv[i]] <- sums[Xv[i]] + Xw[i];
    }
    return sums;
  }
  
  matrix bspline(vector x, int K, int q) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(max(x)-min(x))/(K-q); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,K-q-1)';
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
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Krow; //number of functions in spline base for row biases
  int Kdiag; //number of functions in spline base for diagonal decay
  //biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  int rejoined[S];
  int danglingL[S];
  int danglingR[S];
  //counts : explicit
  int<lower=0> Nclose; //number of close counts modelled explicitly
  int<lower=0> counts_close[Nclose]; //value of the count
  int<lower=0> index_close[2,Nclose]; //indices of rsite pairs
  vector<lower=0>[Nclose] dist_close; //genomic distance between rsites
  //
  int<lower=0> Nfar; //number of far counts modelled explicitly
  int<lower=0> counts_far[Nfar]; //value of the count
  int<lower=0> index_far[2,Nfar]; //indices of rsite pairs
  vector<lower=0>[Nfar] dist_far; //genomic distance between rsites
  //
  int<lower=0> Nup; //number of upstream counts modelled explicitly
  int<lower=0> counts_up[Nup]; //value of the count
  int<lower=0> index_up[2,Nup]; //indices of rsite pairs
  vector<lower=0>[Nup] dist_up; //genomic distance between rsites
  //
  int<lower=0> Ndown; //number of downstream counts modelled explicitly
  int<lower=0> counts_down[Ndown]; //value of the count
  int<lower=0> index_down[2,Ndown]; //indices of rsite pairs
  vector<lower=0>[Ndown] dist_down; //genomic distance between rsites
  //
  //counts : mean field
  int<lower=0> Nl; //number of data points for left side of rsites
  int<lower=0> Nkl_count[Nl]; //value of the count
  int<lower=0> Nkl_cidx[Nl]; //index of that rsite
  int<lower=0> Nkl_N[Nl]; //Nkl(c), i.e. how many contacts have count c in left side of row k
  int<lower=0> Nkl_levels; //number of levels of N
  #
  int<lower=0> Nr; //number of data points for right side of rsites
  int<lower=0> Nkr_count[Nr]; //value of the count
  int<lower=0> Nkr_cidx[Nr]; //index of that rsite
  int<lower=0> Nkr_N[Nr]; //Nkr(c), i.e. how many contacts have count c in right side of row k
  int<lower=0> Nkr_levels; //number of levels of N
  #
  int<lower=0> Nd; //number of data points for mean field on decay
  int<lower=0> Nkd_count[Nd]; //value of the count
  vector[Nd] Nkd_d; //geometric mean distance in bin k
  int<lower=0> Nkd_N[Nd]; //Nkd(c), i.e. how many contacts have count c in decay bin k
  int<lower=0> Nkd_levels; //number of levels of N
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[Krow] prow;
  //diagonal SCAM spline, dense, exact and mean field
  matrix[Nclose+Nfar+Nup+Ndown+Nd,Kdiag] Xdiag;
  row_vector[Kdiag-1] pdiag;
  //indices of changes for N
  int Nkl_Nidx[Nkl_levels+1];
  int Nkr_Nidx[Nkr_levels+1];
  int Nkd_Nidx[Nkd_levels+1];
  
  ////bias spline, sparse (nu and delta have the same design)
  //BEGIN sparse calculation
  //cannot write function that modifies its arguments so we put it here
  //input: vector[Nup+Ndown+Nclose+Nfar] cutsites, int Krow, int splinedegree()
  //output: vector[nnz(S)] Xrow_w, int Xrow_v[nnz(S)], int Xrow_u[S+1]
  {
    real dx; //interval length
    row_vector[Krow] t; //Krownot locations (except last)
    int x_n[Krow-splinedegree()+1]; //cumulative histogram of cutsites values
    int idx_u; //counters for filling of sparse matrix
    int idx_w;
    //
    dx <- 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree()); //maKrowe it slightly larger
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
  prow <- column_sums(Krow, Xrow_w, Xrow_v)';
  prow <- prow / sqrt(dot_self(prow));
  
  //diagonal SCAM spline, dense, exact and mean field model
  {
    //can't do abs() on a vector so be inventive
    vector[Nclose+Nfar+Nup+Ndown] tmp;
    row_vector[Kdiag] tmp2;
    tmp[:Nclose] <- dist_close;
    tmp[(Nclose+1):(Nclose+Nfar)] <- dist_far;
    tmp[(Nclose+Nfar+1):(Nclose+Nfar+Nup)] <- dist_up;
    tmp[(Nclose+Nfar+Nup+1):] <- dist_down;
    Xdiag <- bspline(append_row(log(tmp), log(Nkd_d)), Kdiag, splinedegree());
    //projector for diagonal (SCAM)
    tmp2 <- rep_row_vector(1,Nup+Ndown+Nclose+Nfar+Nd) * Xdiag;
    pdiag <- -tmp2[2:] / (tmp2 * rep_vector(1,Kdiag)); //should it be weighted along MF?
  }
  
  //indices of changes for N
  Nkl_Nidx <- change_points(Nkl_N, Nkl_levels);  
  Nkr_Nidx <- change_points(Nkr_N, Nkr_levels);  
  Nkd_Nidx <- change_points(Nkd_N, Nkd_levels);
  
}
parameters {
  //exposures
  real eC;
  real eRJ;
  real eDE;
  //spline parameters
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  positive_ordered[Kdiag-1] beta_diag;
  //deviances
  real<lower=0> alpha;
  real<lower=0> alpha_mf;
  //length scales
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
  real<lower=0> lambda_diag;
}
transformed parameters {
  //nu
  vector[S] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff; //2nd order difference on beta_nu_aug
  //delta
  vector[S] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff; //2nd order difference on beta_delta_aug
  //diag
  vector[Nclose] log_decay_close;
  vector[Nfar] log_decay_far;
  vector[Nup] log_decay_up;
  vector[Ndown] log_decay_down;
  vector[Nd] log_decay_mf;
  vector[Kdiag-2] beta_diag_diff;
  //means
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  vector[S] log_mean_RJ;
  //
  vector[Nclose] log_mean_cclose;
  vector[Nfar] log_mean_cfar;
  vector[Nup] log_mean_cup;
  vector[Ndown] log_mean_cdown;
  //
  vector[Nl] log_mean_left;
  vector[Nr] log_mean_right;
  vector[Nd] log_mean_decay;
  
  //nu
  {
    vector[Krow] beta_nu_centered;
    vector[Krow] beta_nu_aug;
    beta_nu_aug[1] <- sum(beta_nu);
    beta_nu_aug[2:] <- beta_nu;
    beta_nu_centered <- beta_nu_aug - (beta_nu_aug' * prow) * prow;
    log_nu <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_nu_centered);
    beta_nu_diff <- beta_nu_aug[:(Krow-2)]-2*beta_nu_aug[2:(Krow-1)]+beta_nu_aug[3:];
  }
  //delta
  {
    vector[Krow] beta_delta_centered;
    vector[Krow] beta_delta_aug;
    beta_delta_aug[1] <- sum(beta_delta);
    beta_delta_aug[2:] <- beta_delta;
    beta_delta_centered <- beta_delta_aug - (beta_delta_aug' * prow) * prow;
    log_delta <- csr_matrix_times_vector(S, Krow, Xrow_w, Xrow_v, Xrow_u, beta_delta_centered);
    beta_delta_diff <- beta_delta_aug[:(Krow-2)]-2*beta_delta_aug[2:(Krow-1)]+beta_delta_aug[3:];
  }
  //decay
  {
    vector[Nclose+Nfar+Nup+Ndown+Nd] log_decay;
    vector[Kdiag] beta_diag_centered;
    real epsilon;
    real val;
    epsilon <- -1; //+1 for increasing, -1 for decreasing spline
    val <- epsilon*pdiag*beta_diag;
    beta_diag_centered[1] <- val;
    beta_diag_centered[2:] <- val+epsilon*beta_diag;
    log_decay <- Xdiag * beta_diag_centered;
    log_decay_close <- log_decay[:Nclose];
    log_decay_far <- log_decay[(Nclose+1):(Nclose+Nfar)];
    log_decay_up <- log_decay[(Nclose+Nfar+1):(Nclose+Nfar+Nup)];
    log_decay_down <- log_decay[(Nclose+Nfar+Nup+1):(Nclose+Nfar+Nup+Ndown)];
    log_decay_mf <- log_decay[(Nclose+Nfar+Nup+Ndown+1):];
    beta_diag_diff <- beta_diag_centered[:(Kdiag-2)]-2*beta_diag_centered[2:(Kdiag-1)]+beta_diag_centered[3:];
  }

  //means
  {
    //biases
    log_mean_RJ <- log_nu + eRJ;
    log_mean_DL <- log_nu + eDE + log_delta;
    log_mean_DR <- log_nu + eDE - log_delta;
    //exact counts  
    log_mean_cclose <- eC + log_decay_close + (log_nu - log_delta)[index_close[1]] + (log_nu + log_delta)[index_close[2]];
    log_mean_cfar   <- eC + log_decay_far   + (log_nu + log_delta)[index_far[1]]   + (log_nu - log_delta)[index_far[2]];
    log_mean_cup    <- eC + log_decay_up    + (log_nu + log_delta)[index_up[1]]    + (log_nu + log_delta)[index_up[2]];
    log_mean_cdown  <- eC + log_decay_down  + (log_nu - log_delta)[index_down[1]]  + (log_nu - log_delta)[index_down[2]];
    //mean field counts
    log_mean_left <- eC + (log_nu + log_delta)[Nkl_cidx];
    log_mean_right <- eC + (log_nu - log_delta)[Nkr_cidx];
    log_mean_decay <- eC + log_decay_mf;
  }
}
model {
  //// Exact likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  //counts: Close, Far, Up, Down
  counts_close ~ neg_binomial_2_log(log_mean_cclose, alpha); // Close
  counts_far   ~ neg_binomial_2_log(log_mean_cfar, alpha); // Far
  counts_up    ~ neg_binomial_2_log(log_mean_cup, alpha); // Up
  counts_down  ~ neg_binomial_2_log(log_mean_cdown, alpha); // Down
  
  //// Mean field likelihoods
  // weights must sum to 1, here they are 1/3 1/3 1/3
  //Left
  for (i in 1:Nkl_levels) {
    int begin;
    int end;
    begin <- Nkl_Nidx[i];
    end <- Nkl_Nidx[i+1]-1;
    increment_log_prob(Nkl_N[begin]*neg_binomial_2_log_log(Nkl_count[begin:end], log_mean_left[begin:end], alpha_mf)/3);
  }
  //Right
  for (i in 1:Nkr_levels) {
    int begin;
    int end;
    begin <- Nkr_Nidx[i];
    end <- Nkr_Nidx[i+1]-1;
    increment_log_prob(Nkr_N[begin]*neg_binomial_2_log_log(Nkr_count[begin:end], log_mean_right[begin:end], alpha_mf)/3);
  }
  //Decay
  for (i in 1:Nkd_levels) {
    int begin;
    int end;
    begin <- Nkd_Nidx[i];
    end <- Nkd_Nidx[i+1]-1;
    increment_log_prob(Nkd_N[begin]*neg_binomial_2_log_log(Nkd_count[begin:end], log_mean_decay[begin:end], alpha_mf)/3);
  }
  
  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_nu_diff ~ normal(0, 1./lambda_nu);
  beta_delta_diff ~ normal(0, 1./lambda_delta);
  beta_diag_diff ~ normal(0, 1./lambda_diag);
}
generated quantities {
  real deviance;
  real deviance_null;
  real deviance_proportion_explained;
  #
  real rejoined_deviance;
  real rejoined_deviance_null;
  real rejoined_deviance_proportion_explained;
  #
  real dangling_deviance;
  real dangling_deviance_null;
  real dangling_deviance_proportion_explained;
  #
  real count_deviance_ex;
  real count_deviance_ex_null;
  real count_deviance_ex_proportion_explained;
  #
  real count_deviance_mf;
  real count_deviance_mf_null;
  real count_deviance_mf_proportion_explained;
  #deviances
  count_deviance_mf <- neg_binomial_2_log_deviance(Nkl_count, log_mean_left, alpha, to_vector(Nkl_N)) +
                       neg_binomial_2_log_deviance(Nkr_count, log_mean_right, alpha, to_vector(Nkr_N)) +
                       neg_binomial_2_log_deviance(Nkd_count, log_mean_decay, alpha, to_vector(Nkd_N));
  count_deviance_ex <- neg_binomial_2_log_deviance(counts_close, log_mean_cclose, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_far, log_mean_cfar, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_up, log_mean_cup, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_down, log_mean_cdown, alpha, rep_vector(1,1));
  rejoined_deviance <- neg_binomial_2_log_deviance(rejoined, log_mean_RJ, alpha, rep_vector(1,1));
  dangling_deviance <- neg_binomial_2_log_deviance(danglingL, log_mean_DL, alpha, rep_vector(1,1)) +
                       neg_binomial_2_log_deviance(danglingR, log_mean_DR, alpha, rep_vector(1,1));
  deviance <- rejoined_deviance + dangling_deviance + count_deviance_ex + count_deviance_mf;
  #null deviances
  rejoined_deviance_null <- neg_binomial_2_log_deviance(rejoined, rep_vector(eRJ,1), alpha, rep_vector(1,1));
  dangling_deviance_null <- neg_binomial_2_log_deviance(danglingL, rep_vector(eDE,1), alpha, rep_vector(1,1)) +
                            neg_binomial_2_log_deviance(danglingR, rep_vector(eDE,1), alpha, rep_vector(1,1));
  count_deviance_ex_null <- neg_binomial_2_log_deviance(counts_close, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_far, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_up, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_down, rep_vector(eC,1), alpha, rep_vector(1,1));
  count_deviance_mf_null <- neg_binomial_2_log_deviance(Nkl_count, rep_vector(eC,1), alpha, to_vector(Nkl_N)) +
                            neg_binomial_2_log_deviance(Nkr_count, rep_vector(eC,1), alpha, to_vector(Nkr_N)) +
                            neg_binomial_2_log_deviance(Nkd_count, rep_vector(eC,1), alpha, to_vector(Nkd_N));
  deviance_null <- rejoined_deviance_null + dangling_deviance_null + count_deviance_ex_null + count_deviance_mf_null;
  #proportions explained
  rejoined_deviance_proportion_explained <- 100*(rejoined_deviance_null - rejoined_deviance)/rejoined_deviance_null;
  dangling_deviance_proportion_explained <- 100*(dangling_deviance_null - dangling_deviance)/dangling_deviance_null;
  count_deviance_ex_proportion_explained <- 100*(count_deviance_ex_null - count_deviance_ex)/count_deviance_ex_null;
  count_deviance_mf_proportion_explained <- 100*(count_deviance_mf_null - count_deviance_mf)/count_deviance_mf_null;
  deviance_proportion_explained <- 100*(deviance_null - deviance)/deviance_null;
}
