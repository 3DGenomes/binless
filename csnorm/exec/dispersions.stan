////
// Cut-site normalization model: predict dispersion on binned matrices
////
data {
  //dimensions
  int<lower=1> B;
  //vectorized matrices
  int<lower=0> observed[B]; //summed counts per bin
  vector<lower=0>[B] expected; //posterior mean of negative binomial per bin
  vector<lower=0>[B] ncounts; //number of possible counts in bin
}
parameters {
  real<lower=0> alpha; //dispersion estimate per count
}
transformed parameters {
  vector<lower=0>[B] dispersion;
  dispersion = alpha*ncounts;
}
model {
  observed ~ neg_binomial_2(expected,dispersion);
}
