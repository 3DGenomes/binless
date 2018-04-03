#ifndef FAST_ESTIMATOR_HPP
#define FAST_ESTIMATOR_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "util.hpp" //bin_data_evenly
#include "fast_residuals_pair.hpp"

namespace binless {
namespace fast {

template<class EstimatorImpl>
class Estimator : private EstimatorImpl {
public:
  template<typename FastData, typename Config>
  Estimator(const FastData& data, const Config& conf) : EstimatorImpl(data,conf) {}
  
  //compute group sums of a vector of the size of the input data into the bins formed for the decay calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return EstimatorImpl::get_settings().get_binner()*vec; }
  
  //get log decay along binned distances
  Eigen::VectorXd get_binned_estimate() const {
    return EstimatorImpl::get_estimate() - Eigen::VectorXd::Constant(EstimatorImpl::get_settings().get_nbins(), EstimatorImpl::get_params().get_mean());
  }
  
  //get approximate log decay along distances in original data (same approx as during fitting)
  Eigen::VectorXd get_data_estimate() const {
    return EstimatorImpl::get_settings().get_binner().transpose()*get_binned_estimate();
  }
  
  //provide a way to store and recall the state of the estimator
  Rcpp::List get_state() const { return EstimatorImpl::get_params().get_state(); }
  void set_state(const Rcpp::List& state) { EstimatorImpl::get_params().set_state(state); }
  
  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const std::vector<double>& log_expected, const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z);
  //perform spline fit of summary data and center final estimate
  void update_params() {
    EstimatorImpl::update_params();
    center_estimate(); 
  }
  
private:
  //compute average of estimate in order to center it
  void center_estimate() {
    EstimatorImpl::get_params().set_mean(EstimatorImpl::get_settings().get_nobs().dot(EstimatorImpl::get_settings().get_X() * EstimatorImpl::get_params().get_beta())
                                           /EstimatorImpl::get_settings().get_nobs().sum());
  }
  
};

#include "fast_estimator.ipp"

}
}

#endif

