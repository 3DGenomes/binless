#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_bias_Mean.hpp"
#include "fast_residuals.hpp"
#include "spline.hpp"
#include "fast_expected.hpp"

namespace binless {
namespace fast {


//here, log bias is log ( sum_i observed / sum_i expected ) with i summed over counter diagonals
void BiasMeanEstimator::set_poisson_lsq_summary(const std::vector<double>& log_expected_std, const FastSignalData& data, double pseudocount) {
  //compute observed data
  auto observed_std = data.get_observed();
  auto nobs_std = data.get_nobs();
  Eigen::VectorXd observed = Eigen::VectorXd::Zero(observed_std.size());
  Eigen::VectorXd nobs = Eigen::VectorXd::Zero(nobs_std.size());
  for (unsigned i=0; i<observed.rows(); ++i) observed(i) = observed_std[i]; // cast to double
  for (unsigned i=0; i<observed.rows(); ++i) nobs(i) = nobs_std[i]; // cast to double
  Eigen::VectorXd sum_obs = summarize((observed.array()*nobs.array()).matrix());
  //compute expected data
  const Eigen::Map<const Eigen::VectorXd> log_expected(log_expected_std.data(),log_expected_std.size());
  Eigen::VectorXd sum_exp = summarize((log_expected.array().exp()*nobs.array()).matrix());
  //compute log_bias
  Eigen::VectorXd log_bias = (pseudocount + sum_obs.array()/sum_exp.array()).log().matrix();
  //center log_bias
  double avg = settings_.get_nobs().dot(log_bias)/settings_.get_nobs().sum();
  log_bias = (log_bias.array() - avg).matrix();
  //compute weight
  Eigen::VectorXd weight = (log_bias.array().exp()*settings_.get_nobs().array()).matrix();
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance etahat weight sum_obs sum_exp\n";
  Rcpp::Rcout << (Eigen::MatrixXd(log_bias.rows(),5) << settings_.get_log_distance().array().exp().matrix(),
                  log_bias, weight, sum_obs, sum_exp).finished();*/
  //report back
  summary_.set_etahat(log_bias);
  summary_.set_weight(weight);
}

void BiasMeanEstimator::update_summary(const ResidualsPair& z) {
  //compute weight
  const Eigen::Map<const Eigen::VectorXd> weights(z.weights.data(),z.weights.size());
  Eigen::VectorXd weight_sum = summarize(weights);
  //compute etahat
  const Eigen::Map<const Eigen::VectorXd> residuals(z.residuals.data(),z.residuals.size());
  Eigen::VectorXd etahat = summarize(residuals.array() * weights.array()).matrix();
  Eigen::VectorXd log_bias = get_binned_estimate();
  etahat = (etahat.array() / weight_sum.array()).matrix() + log_bias;
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance etahat weight\n";
  Rcpp::Rcout << (Eigen::MatrixXd(etahat.rows(),3) << settings_.get_log_distance().array().exp().matrix(), etahat,
                  weight_sum).finished();*/
  summary_.set_etahat(etahat);
  summary_.set_weight(weight_sum);
}

void BiasMeanEstimator::update_params() {
  //center log_bias (no smoothing)
  double avg = settings_.get_nobs().dot(summary_.get_etahat())/settings_.get_nobs().sum();
  Eigen::VectorXd log_biases = (summary_.get_etahat().array() - avg).matrix();
  //cap estimates at 3SD from the mean
  double stdev = std::sqrt(log_biases.squaredNorm()/log_biases.rows());
  log_biases = log_biases.array().min(3*stdev).max(-3*stdev).matrix();
  
  set_estimate(log_biases);
  //center log bias
  //center_estimate();
}

}
}

