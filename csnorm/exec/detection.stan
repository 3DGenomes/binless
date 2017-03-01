////
// Cut-site normalization: interaction detection through model comparison
////
data {
  int<lower=1> G; //number of groups
  int<lower=1> N; // Number of counts in this bin
  int<lower=1,upper=N+1> cbegin[G+1]; //cbegin[i]=j: group i starts at j
  int<lower=0> observed[N]; // Value of each count
  vector[N] log_expected; // Model log mean for each count
  real<lower=0> alpha; //dispersion
  real<lower=0> sigma; //for normal prior
}
transformed data {
  vector[G] bd;
  for (g in 1:G) bd[g] = sum(observed[cbegin[g]:(cbegin[g+1]-1)]);
}
parameters {
  real log_s[G]; //log(signal)
}
model {
  for (g in 1:G) {
    target += neg_binomial_2_log_lpmf(observed[cbegin[g]:(cbegin[g+1]-1)] | log_s[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
  }
  target += normal_lpdf(log_s | 0, sigma);
}
generated quantities {
  vector[G] lpdfs;
  vector[G] lpdf0;
  vector[G] binned;
  for (g in 1:G) {
    lpdfs[g] = neg_binomial_2_log_lpmf(observed[cbegin[g]:(cbegin[g+1]-1)] | log_s[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha)
              + normal_lpdf(log_s[g] | 0, sigma);
    lpdf0[g] = neg_binomial_2_log_lpmf(observed[cbegin[g]:(cbegin[g+1]-1)] | log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
  }
  binned=bd;
}
