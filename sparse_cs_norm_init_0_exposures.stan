////
// Cut-site normalization model: guess eDE eC and eRJ
////
data {
  //biases
  int<lower=1> S; //number of cut sites
  int rejoined[S];
  int danglingL[S];
  int danglingR[S];
  //counts
  int<lower=1> N; //number of data points
  int<lower=0> counts[4,N]; //raw counts: Close, Far, Up, Down
}
parameters {
  real eC;  //exposure for counts
  real eRJ; //exposure for rejoined ends
  real eDE; //exposure for dangling ends
  real<lower=0> alpha;
}
model {
  //// likelihoods
  //biases
  rejoined ~ neg_binomial_2_log(eRJ, alpha);
  danglingL ~ neg_binomial_2_log(eDE, alpha);
  danglingR ~ neg_binomial_2_log(eDE, alpha);
  //counts: Close, Far, Up, Down
  counts[1] ~ neg_binomial_2_log(eC, alpha); // Close
  counts[2] ~ neg_binomial_2_log(eC, alpha); // Far
  counts[3] ~ neg_binomial_2_log(eC, alpha); // Up
  counts[4] ~ neg_binomial_2_log(eC, alpha); // Down
}