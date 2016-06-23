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
}
////////////////////////////////////////////////////////////////
data {
  //genomic biases
  int<lower=4> Krow; //number of functions in spline base for row biases
  int<lower=1> S; //number of cut sites
  vector[S] cutsites; //cut site locations
  int rejoined[S];
  int danglingL[S];
  int danglingR[S];
  //decay bias
  int<lower=4> Kdiag; //number of functions in spline base for diagonal decay
  real<lower=0> dmin; //distance bounds for spline
  real<lower=dmin> dmax;
  //counts : explicit
  int<lower=0> N; //number of counts
  int<lower=1,upper=S> cidx[2,N]; //indices of rsite pairs
  vector<lower=dmin,upper=dmax>[N] dist; //genomic distance between rsites
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;
  //diagonal SCAM spline, dense
  matrix[N,Kdiag] Xdiag;
  row_vector[Kdiag] pdiag;
  vector[N] diag_weights;
  
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
  {
    row_weights <- rep_vector(3, S); //dangling L/R + rejoined
    for (i in 1:N) {
      row_weights[cidx[1,i]] <- row_weights[cidx[1,i]] + 1;
      row_weights[cidx[2,i]] <- row_weights[cidx[2,i]] + 1;
    }
    //
    row_weights <- row_weights/mean(row_weights);
    prow <- vector_times_csr_matrix(Krow, row_weights, Xrow_w, Xrow_v, Xrow_u);
    prow <- prow / (prow * rep_vector(1,Krow));
  }
  
  //diagonal SCAM spline, dense, exact and mean field model
  {
    Xdiag <- bspline(log(to_vector(dist)), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    diag_weights <- rep_vector(1, N);
    pdiag <- diag_weights' * Xdiag;
    pdiag <- pdiag / (pdiag * rep_vector(1,Kdiag));
  }
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
  //dispersion
  real<lower=0> alpha;
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
  vector[N] log_decay;
  vector[Kdiag] beta_diag_centered;
  vector[Kdiag-2] beta_diag_diff;
  //means
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  vector[S] log_mean_RJ;
  //
  vector[N] log_mean_cclose;
  vector[N] log_mean_cfar;
  vector[N] log_mean_cup;
  vector[N] log_mean_cdown;
  
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
    //biases
    log_mean_RJ <- log_nu + eRJ;
    log_mean_DL <- log_nu + eDE + log_delta;
    log_mean_DR <- log_nu + eDE - log_delta;
    //exact counts  
    log_mean_cclose <- eC + log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cfar   <- eC + log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
    log_mean_cup    <- eC + log_decay + (log_nu + log_delta)[cidx[1]] + (log_nu + log_delta)[cidx[2]];
    log_mean_cdown  <- eC + log_decay + (log_nu - log_delta)[cidx[1]] + (log_nu - log_delta)[cidx[2]];
  }
}
model {
  //// likelihoods
  //biases
  rejoined  ~ neg_binomial_2_log(log_mean_RJ, alpha);
  danglingL ~ neg_binomial_2_log(log_mean_DL, alpha);
  danglingR ~ neg_binomial_2_log(log_mean_DR, alpha);
  
  //counts: Close, Far, Up, Down
  counts_close ~ neg_binomial_2_log(log_mean_cclose, alpha); // Close
  counts_far   ~ neg_binomial_2_log(log_mean_cfar, alpha); // Far
  counts_up    ~ neg_binomial_2_log(log_mean_cup, alpha); // Up
  counts_down  ~ neg_binomial_2_log(log_mean_cdown, alpha); // Down
  
  //// Priors
  //P-spline prior on the differences (K-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_nu_diff ~ normal(0, 1./lambda_nu);
  beta_delta_diff ~ normal(0, 1./lambda_delta);
  beta_diag_diff ~ normal(0, 1./lambda_diag);
}
generated quantities {
  vector[S] rw;
  vector[N] dw;
  #
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
  real count_deviance;
  real count_deviance_null;
  real count_deviance_proportion_explained;
  #weights
  rw <- row_weights;
  dw <- diag_weights;
  #deviances
  count_deviance <- neg_binomial_2_log_deviance(counts_close, log_mean_cclose, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_far, log_mean_cfar, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_up, log_mean_cup, alpha, rep_vector(1,1)) +
              neg_binomial_2_log_deviance(counts_down, log_mean_cdown, alpha, rep_vector(1,1));
  rejoined_deviance <- neg_binomial_2_log_deviance(rejoined, log_mean_RJ, alpha, rep_vector(1,1));
  dangling_deviance <- neg_binomial_2_log_deviance(danglingL, log_mean_DL, alpha, rep_vector(1,1)) +
                       neg_binomial_2_log_deviance(danglingR, log_mean_DR, alpha, rep_vector(1,1));
  deviance <- rejoined_deviance + dangling_deviance + count_deviance;
  #null deviances
  rejoined_deviance_null <- neg_binomial_2_log_deviance(rejoined, rep_vector(eRJ,1), alpha, rep_vector(1,1));
  dangling_deviance_null <- neg_binomial_2_log_deviance(danglingL, rep_vector(eDE,1), alpha, rep_vector(1,1)) +
                            neg_binomial_2_log_deviance(danglingR, rep_vector(eDE,1), alpha, rep_vector(1,1));
  count_deviance_null <- neg_binomial_2_log_deviance(counts_close, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_far, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_up, rep_vector(eC,1), alpha, rep_vector(1,1)) +
                   neg_binomial_2_log_deviance(counts_down, rep_vector(eC,1), alpha, rep_vector(1,1));
  deviance_null <- rejoined_deviance_null + dangling_deviance_null + count_deviance_null;
  #proportions explained
  rejoined_deviance_proportion_explained <- 100*(rejoined_deviance_null - rejoined_deviance)/rejoined_deviance_null;
  dangling_deviance_proportion_explained <- 100*(dangling_deviance_null - dangling_deviance)/dangling_deviance_null;
  count_deviance_proportion_explained <- 100*(count_deviance_null - count_deviance)/count_deviance_null;
  deviance_proportion_explained <- 100*(deviance_null - deviance)/deviance_null;
}
