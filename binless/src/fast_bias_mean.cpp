#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_bias_mean.hpp"
#include "fast_residuals_pair.hpp"
#include "spline.hpp"
#include "fast_expected.hpp"

namespace binless {
namespace fast {

void BiasMeanImpl::update_params() {
  //center log_bias (no smoothing)
  double avg = settings_.get_nobs().dot(summary_.get_phihat())/settings_.get_nobs().sum();
  Eigen::VectorXd log_biases = (summary_.get_phihat().array() - avg).matrix();
  //cap estimates at 3SD from the mean
  double stdev = std::sqrt(log_biases.squaredNorm()/log_biases.rows());
  log_biases = log_biases.array().min(3*stdev).max(-3*stdev).matrix();
  
  set_estimate(log_biases);
}

}
}

