////
// Cut-site normalization: binning matrices
////
data {
  int<lower=1> G; //number of groups
  int<lower=1> N; // Number of counts in this bin
  int<lower=1,upper=N+1> cbegin[G+1]; //cbegin[i]=j: group i starts at j
  int<lower=0> count[N]; // Value of each count
  vector[N] log_expected; // Model log mean for each count
  vector[N] log_decay; // Model log decay for each count
  real<lower=0> alpha; //dispersion
}
parameters {
  real log_s[G]; //log(signal)
  real<lower=0> gamma;
}
model {
  for (g in 1:G) {
    target += neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                        log_s[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
  }
  log_s ~ double_exponential(0, gamma);
}
