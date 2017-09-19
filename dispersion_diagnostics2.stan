data {
  int<lower=1> ncounts;
  int<lower=0> count[ncounts];
  vector[ncounts] log_mean;
  int<lower=1> ntypes;
  int<lower=1,upper=ntypes> type[ncounts];
  vector[ncounts] indep;
}
parameters {
  real exposure[ntypes];
  real a;
  real<lower=0> b;
}
transformed parameters {}
model {
  vector[ncounts] exp_vec;
  vector[ncounts] disp_vec;
  for (i in 1:ncounts) {
    exp_vec[i] = exposure[type[i]];
  }
  disp_vec = b*exp(a*indep);
  count  ~ neg_binomial_2_log(log_mean + exp_vec, disp_vec);
}
