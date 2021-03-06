#ifndef FAST_FITTER_HPP
#define FAST_FITTER_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "fast_fitter_gam.hpp"
#include "fast_fitter_mean.hpp"

namespace binless {
namespace fast {

template<typename Leg, typename Method>
class Fitter : private FitterImpl<Leg,Method> {
public:
  typedef FitterImpl<Leg,Method> fitterImpl_t;
  
  Fitter(const SummarizerSettings& settings, const Config<Leg,Method>& conf) : fitterImpl_t(settings,conf) {}

  //provide a way to store and recall the state of the Fitter
  Rcpp::List get_state() const { return fitterImpl_t::get_params().get_state(); }
  void set_state(const Rcpp::List& state) { fitterImpl_t::get_params().set_state(state); }
  
  //perform fit of summary data
  void update_params(const Eigen::VectorXd& phihat, const Eigen::VectorXd& weight) {
    fitterImpl_t::update_params(phihat,weight);
  }

  //get estimate along binned support
  Eigen::VectorXd get_binned_estimate() const {
    return fitterImpl_t::get_estimate() - Eigen::VectorXd::Constant(fitterImpl_t::get_settings().get_nbins(), fitterImpl_t::get_params().get_mean());
  }
  
};

}
}

#endif

