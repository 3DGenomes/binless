#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_bias.hpp"
#include "fast_residuals.hpp"
#include "spline.hpp"
#include "gam.hpp"
#include "fast_expected.hpp"

namespace binless {
namespace fast {


//here, log bias is log ( sum_i observed / sum_i expected ) with i summed over counter diagonals
void BiasEstimator::set_poisson_lsq_summary(const FastSignalData& data, double pseudocount) {
  //compute observed data
  auto observed_std = data.get_observed();
  auto nobs_std = data.get_nobs();
  Eigen::VectorXd observed = Eigen::VectorXd::Zero(observed_std.size());
  Eigen::VectorXd nobs = Eigen::VectorXd::Zero(nobs_std.size());
  for (unsigned i=0; i<observed.rows(); ++i) observed(i) = observed_std[i]; // cast to double
  for (unsigned i=0; i<observed.rows(); ++i) nobs(i) = nobs_std[i]; // cast to double
  Eigen::VectorXd sum_obs = summarize((observed.array()*nobs.array()).matrix());
  //compute expected data
  auto log_expected_std = get_log_expected(data, *this);
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
  summary_.etahat = log_bias;
  summary_.weight = weight;
}

void BiasEstimator::update_summary(const ResidualsPair& z) {
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
  summary_.etahat = etahat;
  summary_.weight = weight_sum;
}

void BiasEstimator::update_params() {
  //extract data
  const Eigen::VectorXd y(summary_.etahat);
  const Eigen::VectorXd Sm1 = summary_.weight.array().sqrt().matrix();
  
  gam_.optimize(y,Sm1,settings_.get_max_iter(),settings_.get_tol_val());
  //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
  set_lambda(gam_.get_lambda());
  set_beta(gam_.get_beta());
  //center log bias
  center_estimate();
  /*Rcpp::Rcout << "spline_log_bias_fit\n";
  Rcpp::Rcout << "distance etahat weight nobs log_bias\n";
  Rcpp::Rcout << (Eigen::MatrixXd(settings_.get_nbins(),5) << settings_.get_log_distance().array().exp().matrix(),
                  summary_.etahat, summary_.weight, settings_.get_nobs(), get_binned_estimate()).finished();*/
}

}
}

