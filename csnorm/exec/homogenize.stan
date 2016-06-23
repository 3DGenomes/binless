////
// Cut-site normalization model: retrieve eRJ eDE log_nu and log_delta
//                               at constant stiffness
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int<lower=0> Krow; //number of functions in spline base for row biases
  //biases
  int<lower=1> S;
  vector[S] cutsites;
  vector[S] log_mean_RJ;
  vector[S] log_mean_DL;
  vector[S] log_mean_DR;
  //stiffnesses
  real<lower=0> lambda_nu;
  real<lower=0> lambda_delta;
}
transformed data {
  //bias spline, sparse (nu and delta have the same design)
  vector[nnz(S)] Xrow_w;
  int Xrow_v[nnz(S)];
  int Xrow_u[S+1];
  vector[S] row_weights;
  row_vector[Krow] prow;
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
  //build projector on centering constraint
  prow <- column_sums(Krow, Xrow_w, Xrow_v);
  prow <- prow / (prow * rep_vector(1,Krow));
}
parameters {
  //exposures
  real eRJ;
  real eDE;
  //spline coefs
  vector[Krow-1] beta_nu;
  vector[Krow-1] beta_delta;
  //SD
  real<lower=0,upper=1> sigma;
}
transformed parameters {
  //nu
  vector[S] log_nu; // log(nu)
  vector[Krow-2] beta_nu_diff; //2nd order difference on beta_nu_aug
  //delta
  vector[S] log_delta; // log(delta)
  vector[Krow-2] beta_delta_diff; //2nd order difference on beta_delta_aug
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
}
model {
  //Likelihoods
  log_mean_RJ ~ normal(eRJ + log_nu, sigma);
  log_mean_DL ~ normal(eDE  + log_nu + log_delta, sigma);
  log_mean_DR ~ normal(eDE  + log_nu - log_delta, sigma);
  //P-spline prior on the 2nd order differences (Krow-2 params)
  //warning on jacobian can be ignored
  //see GAM, Wood (2006), section 4.8.2 (p.187)
  beta_nu_diff ~ normal(0, 1/lambda_nu);
  beta_delta_diff ~ normal(0, 1/lambda_delta);
}
