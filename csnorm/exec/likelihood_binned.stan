////
// Cut-site normalization: evaluate signal likelihood at a given value
////
data {
  int<lower=1> G; //number of groups
  int<lower=1> N; // Number of counts in this bin
  int<lower=1,upper=N+1> cbegin[G+1]; //cbegin[i]=j: group i starts at j
  int<lower=0> count[N]; // Value of each count
  vector[N] log_expected; // Model log mean for each count
  real<lower=0> alpha; //dispersion
  real log_s[G]; //log(signal)
  int<lower=1> S; //number of intermediate signal positions
}
parameters {}
model {}
generated quantities {
  vector[G] lpdfref;
  matrix[G,S] lpdfs;
  vector[G] lpdf0;

  for (g in 1:G) {
    lpdfref[g] = neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                        log_s[g]+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
    for (s in 1:S) {
      lpdfs[g,s] = neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                        (log_s[g]*(s-1))/(S-1)+log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
    }
    lpdf0[g] = neg_binomial_2_log_lpmf(count[cbegin[g]:(cbegin[g+1]-1)] |
                      log_expected[cbegin[g]:(cbegin[g+1]-1)], alpha);
  }
}
