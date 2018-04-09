#ifndef FAST_SUMMARIZER_HPP
#define FAST_SUMMARIZER_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "util.hpp" //bin_data_evenly
#include "fast_residuals_pair.hpp"
#include "FastData.hpp"
#include "macros.hpp"
#include "Traits.hpp"

namespace binless {
namespace fast {

// class that holds summary statistics (aka IRLS weights)
struct Summary {
  Summary() : phihat_(Eigen::VectorXd()), weight_(Eigen::VectorXd()) {}
  BINLESS_FORBID_COPY(Summary);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, phihat);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, weight);
};

template<typename Leg, typename Method>
class SummarizerImpl {
public:
  template<typename FastData>
  SummarizerImpl(const FastData& data, const Config<Leg,Method>& conf) : 
    settings_(data, conf), summary_() {}
  
  BINLESS_GET_CONSTREF_DECL(SummarizerSettings<Leg COMMA() Method>, settings);
  BINLESS_GET_REF_DECL(Summary, summary);
  
};

// Estimator is a policy class that takes data and some configuration info and performs a complete IRLS step on it
template<class Leg, typename Method>
class Summarizer : public SummarizerImpl<Leg,Method> {
public:
  typedef SummarizerImpl<Leg,Method> summarizerImpl_t;
  
  template<typename FastData, typename Config>
  Summarizer(const FastData& data, const Config& conf) : summarizerImpl_t(data,conf) {}
  
  //compute group sums of a vector of the size of the input data into the bins formed for the decay calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return summarizerImpl_t::get_settings().get_binner()*vec; }
  
  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const std::vector<double>& log_expected, const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z, const Eigen::VectorXd& estimate);
  
};

#include "fast_summarizer.ipp"

}
}

#endif

