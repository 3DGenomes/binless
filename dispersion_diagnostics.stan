data {
  int<lower=1> ncounts;
  int<lower=0> count[ncounts];
  vector[ncounts] log_mean;
  int<lower=1> ntypes;
  int<lower=1,upper=ntypes> type[ncounts];
  int<lower=1> ngroups;
  int<lower=1,upper=ngroups> dispgroup[ncounts];
}
parameters {
  real exposure[ntypes];
  real<lower=0> dispersion[ngroups];
}
model {
  vector[ncounts] exp_vec;
  vector[ncounts] disp_vec;
  for (i in 1:ncounts) {
    exp_vec[i] = exposure[type[i]];
    disp_vec[i] = dispersion[dispgroup[i]];
  }
  count  ~ neg_binomial_2_log(log_mean + exp_vec, disp_vec);
}
