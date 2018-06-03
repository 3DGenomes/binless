#ifndef FAST_ESTIMATOR_HPP
#define FAST_ESTIMATOR_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "util.hpp" //bin_data_evenly
#include "fast_residuals_pair.hpp"
#include "FastData.hpp"
#include "macros.hpp"
#include "fast_summarizer.hpp"
#include "fast_fitter.hpp"

namespace binless {
namespace fast {

template<typename Leg, typename Method>
class Estimator : public Summarizer<Leg,Method>, public Fitter<Leg,Method> {
public:
  
  typedef Summarizer<Leg,Method> summarizer_t;
  typedef Fitter<Leg,Method> fitter_t;
  
  template<typename FastData>
  Estimator(const FastData& data, const Config<Leg,Method>& conf) : summarizer_t(data,conf), fitter_t(summarizer_t::get_settings(),conf) {}
  
  void update_summary(const ResidualsPair& z) {
    Eigen::VectorXd estimate = fitter_t::get_binned_estimate();
    summarizer_t::update_summary(z,estimate);
  }
  
  void update_params() {
    Eigen::VectorXd phihat = summarizer_t::get_summary().get_phihat();
    Eigen::VectorXd weight = summarizer_t::get_summary().get_weight();
    fitter_t::update_params(phihat,weight);
  }

  //get approximate estimate along support in original data (same approx as during fitting)
  Eigen::VectorXd get_data_estimate() const {
    return summarizer_t::get_settings().get_binner().transpose()*fitter_t::get_binned_estimate();
  }

};

}
}

#endif

