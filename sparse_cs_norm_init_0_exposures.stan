////
// Cut-site normalization model: guess eDE eC and eRJ
////
functions {
  real neg_binomial_2_log_deviance(int[] y, vector log_mu, real alpha) {
    vector[size(y)] y_vec;
    vector[size(y)] y_mod;
    vector[rows(log_mu)] mu;
    if (rows(log_mu) != size(y)) reject("sizes of y (",size(y),") and log_mu (", rows(log_mu),
                                        ") must match")
    y_vec <- to_vector(y);
    for (i in 1:size(y)) if (y[i]>0) {y_mod[i] <- y[i];} else {y_mod[i] <-1;}
    mu <- exp(log_mu);
    return 2*sum( (y_vec+alpha) .* log( (mu+alpha) ./ (y_vec+alpha) )
                  + y_vec .* log(y_mod ./ mu));
  }
}
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
generated quantities {
  real deviance;
  deviance <- neg_binomial_2_log_deviance(rejoined, rep_vector(eRJ, S), alpha)+
              neg_binomial_2_log_deviance(danglingL, rep_vector(eDE, S), alpha)+
              neg_binomial_2_log_deviance(danglingR, rep_vector(eDE, S), alpha)+
              neg_binomial_2_log_deviance(counts[1], rep_vector(eC, N), alpha)+
              neg_binomial_2_log_deviance(counts[2], rep_vector(eC, N), alpha)+
              neg_binomial_2_log_deviance(counts[3], rep_vector(eC, N), alpha)+
              neg_binomial_2_log_deviance(counts[4], rep_vector(eC, N), alpha);
}
