data {
  int<lower=1> N; // Number of counts in this bin
  int<lower=0> observed[N]; // Value of each count
  vector[N] log_expected; // Model log mean for each count
  real<lower=0> alpha; //dispersion
}
parameters {
  real log_s; //log(signal)
}
model {
  observed ~ neg_binomial_2_log(log_s+log_expected, alpha);
}
