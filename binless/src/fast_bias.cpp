#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_bias.hpp"

namespace binless {
namespace fast {

void BiasGAMImpl::update_params() {
  //extract data
  const Eigen::VectorXd y(summary_.get_phihat());
  const Eigen::VectorXd Sm1 = summary_.get_weight().array().sqrt().matrix();
  gam_.optimize(y,Sm1,settings_.get_max_iter(),settings_.get_tol_val());
  //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
  set_lambda(gam_.get_lambda());
  set_beta(gam_.get_beta());
}

}
}

