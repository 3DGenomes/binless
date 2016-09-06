////
// Cut-site normalization model: predict dispersion on binned matrices
////
data {
  //dimensions
  int<lower=1> B;
  //vectorized matrices
  int<lower=0> observed[B]; //summed counts per bin
  vector<lower=0>[B] expected; //posterior mean of negative binomial per bin
}
parameters {
  real<lower=0> alpha; //dispersion estimate per count
}
model {
  observed ~ neg_binomial_2(expected,alpha);
}
generated quantities {
  vector<lower=0>[B] dispersion;
  dispersion = rep_vector(alpha,B);
}
