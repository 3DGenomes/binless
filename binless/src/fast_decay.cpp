#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_decay.hpp"
#include "fast_residuals.hpp"
#include "spline.hpp"
#include "gam.hpp"
#include "fast_expected.hpp"

namespace binless {
namespace fast {


//here, log decay is log ( sum_i observed / sum_i expected ) with i summed over counter diagonals
void Decay::set_poisson_lsq_summary(const FastSignalData& data, double pseudocount) {
  //compute observed data
  auto observed_std = data.get_observed();
  Eigen::VectorXd observed = Eigen::VectorXd::Zero(observed_std.size());
  for (unsigned i=0; i<observed.rows(); ++i) observed(i) = observed_std[i]; // cast to double
  Eigen::VectorXd sum_obs = summarize(observed);
  //compute expected data
  auto log_expected_std = get_log_expected(data, *this);
  const Eigen::Map<const Eigen::VectorXd> log_expected(log_expected_std.data(),log_expected_std.size());
  Eigen::VectorXd sum_exp = summarize(log_expected.array().exp().matrix());
  //compute log_decay
  Eigen::VectorXd log_decay = ((sum_obs.array() + pseudocount)/sum_exp.array()).log().matrix();
  //center log_decay
  double avg = log_decay.sum()/log_decay.rows();
  log_decay = (log_decay.array() - avg).matrix();
  //compute weight
  Eigen::VectorXd weight = log_decay.array().exp().matrix();
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance kappahat weight sum_obs sum_exp\n";
  Rcpp::Rcout << (Eigen::MatrixXd(log_decay.rows(),5) << schedule_.get_log_distance().array().exp().matrix(),
                  log_decay, weight, sum_obs, sum_exp).finished();*/
  //report back
  summary_.kappahat = log_decay;
  summary_.weight = weight;
}

void Decay::update_summary(const FastSignalData& data) {
  //get residuals
  ResidualsPair z = get_poisson_residuals(data, *this);
  //compute weight
  const Eigen::Map<const Eigen::VectorXd> weights(z.weights.data(),z.weights.size());
  Eigen::VectorXd weight_sum = summarize(weights);
  //compute kappahat
  const Eigen::Map<const Eigen::VectorXd> residuals(z.residuals.data(),z.residuals.size());
  Eigen::VectorXd kappahat = summarize(residuals.array() * weights.array()).matrix();
  Eigen::VectorXd log_decay = get_binned_log_decay();
  kappahat = (kappahat.array() / weight_sum.array()).matrix() + log_decay;
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance kappahat weight\n";
  Rcpp::Rcout << (Eigen::MatrixXd(kappahat.rows(),3) << schedule_.get_log_distance().array().exp().matrix(), kappahat,
                  weight_sum).finished();*/
  summary_.kappahat = kappahat;
  summary_.weight = weight_sum;
}

void Decay::update_params() {
  //extract data
  const Eigen::VectorXd y(summary_.kappahat);
  const Eigen::VectorXd Sm1 = summary_.weight.array().sqrt().matrix();
  
  gam_.optimize(y,Sm1,schedule_.get_max_iter(),schedule_.get_tol_val());
  //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
  set_lambda_diag(gam_.get_lambda());
  set_beta_diag(gam_.get_beta());
  //center log decay
  center_log_decay();
  /*Rcpp::Rcout << "spline_log_decay_fit\n";
  Rcpp::Rcout << "distance kappahat weight ncounts log_decay\n";
  Rcpp::Rcout << (Eigen::MatrixXd(schedule_.get_nbins(),5) << schedule.get_log_distance().array().exp().matrix(),
                  summary_.kappahat, summary_.weight, schedule_.get_ncounts(), get_binned_log_decay()).finished();*/
}

}
}

