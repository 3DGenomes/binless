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

template<typename SummarizerImpl, typename FitterImpl>
class Estimator : public Summarizer<SummarizerImpl>, public Fitter<FitterImpl> {
public:
  
  typedef Summarizer<SummarizerImpl> summarizer_t;
  typedef Fitter<FitterImpl> fitter_t;
  
  template<typename FastData, typename Config>
  Estimator(const FastData& data, const Config& conf) : summarizer_t(data,conf), fitter_t(summarizer_t::get_settings(),data,conf) {}
  
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

