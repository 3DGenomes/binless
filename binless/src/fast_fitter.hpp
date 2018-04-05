#ifndef FAST_FITTER_HPP
#define FAST_FITTER_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "fast_fitter_gam.hpp"
#include "fast_fitter_mean.hpp"

namespace binless {
namespace fast {

template<class FitterImpl>
class Fitter : private FitterImpl {
public:
  template<typename SummarizerSettings, typename FastData, typename Config>
  Fitter(const SummarizerSettings& settings, const FastData& data, const Config& conf) : FitterImpl(settings,data,conf) {}

  //provide a way to store and recall the state of the Fitter
  Rcpp::List get_state() const { return FitterImpl::get_params().get_state(); }
  void set_state(const Rcpp::List& state) { FitterImpl::get_params().set_state(state); }
  
  //perform fit of summary data and center final estimate
  void update_params(const Eigen::VectorXd& phihat, const Eigen::VectorXd& weight) {
    FitterImpl::update_params(phihat,weight);
    if (FitterImpl::get_settings().is_centered()) center_estimate(); 
  }

  //get estimate along binned support
  Eigen::VectorXd get_binned_estimate() const {
    return FitterImpl::get_estimate() - Eigen::VectorXd::Constant(FitterImpl::get_settings().get_nbins(), FitterImpl::get_params().get_mean());
  }
  
private:
  //compute average of estimate in order to center it
  void center_estimate() {
    FitterImpl::get_params().set_mean(FitterImpl::get_settings().get_nobs().dot(FitterImpl::get_estimate())
                                           /FitterImpl::get_settings().get_nobs().sum());
  }
};

}
}

#endif

