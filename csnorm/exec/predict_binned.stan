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
  real<lower=0> gamma; //for double_exponential prior
}
parameters {
  real log_s[G]; //log(signal)
  real log_r[G]; //for ice-like matrix
}
model {
  for (g in 1:G) {
    target += neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                        log_s[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
    target += neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                        log_r[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)]-log_decay[cbegin[g]:(cbegin[g+1]-1)], alpha);
  }
  log_s ~ double_exponential(0, gamma);
  log_r ~ double_exponential(0, gamma);
}
generated quantities {
  vector[G] lpdfr;
  vector[G] lpdfs;
  vector[G] lpdf0;
  int ncounts[G];
  int observed[G];
  vector[G] expected;
  vector[G] expected_sd;
  vector[G] decaymat;
  
  for (g in 1:G) {
    lpdfr[g] = neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                      log_r[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)]-log_decay[cbegin[g]:(cbegin[g+1]-1)], alpha)
              + double_exponential_lpdf(log_r[g] | 0, gamma);
    lpdfs[g] = neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                      log_s[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha) + double_exponential_lpdf(log_s[g] | 0, gamma);
    lpdf0[g] = neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                      log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
    ncounts[g] = cbegin[g+1]-cbegin[g];
    observed[g]=0;
    expected[g]=0;
    expected_sd[g]=0;
    decaymat[g]=0;
    for (i in cbegin[g]:(cbegin[g+1]-1)) {
      real tmp;
      observed[g] = observed[g] + count[i];
      tmp = exp(log_expected[i]);
      expected[g] = expected[g] + tmp;
      expected_sd[g] = expected_sd[g] + tmp + tmp*(tmp/alpha);
      decaymat[g] = decaymat[g] + exp(log_expected[i] - log_decay[i]);
    }
    expected_sd[g] = sqrt(expected_sd[g]);
  }
}
