data {
  int<lower=1> N;
  int<lower=0> observed[N];
}
parameters {
  real expected;
  real<lower=0> alpha;
}
model {
  observed ~ neg_binomial_2(expected,alpha);
}
