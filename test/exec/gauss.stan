data {
  int<lower=1> N;
  vector[N] observed;
}
parameters {
  real expected;
  real<lower=0> alpha;
}
model {
  observed ~ normal(expected,alpha);
}
