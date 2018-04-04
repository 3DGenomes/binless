#ifndef FAST_ESTIMATOR_HPP
#define FAST_ESTIMATOR_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "util.hpp" //bin_data_evenly
#include "fast_residuals_pair.hpp"
#include "FastData.hpp"
#include "macros.hpp"

namespace binless {
namespace fast {

// class that holds summary statistics (aka IRLS weights)
struct Summary {
  Summary() : phihat_(Eigen::VectorXd()), weight_(Eigen::VectorXd()) {}
  BINLESS_FORBID_COPY(Summary);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, phihat);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, weight);
};

// class that holds parameters that were output by a GAM fit
struct GAMParams {
  template<typename Settings>
  GAMParams(const Settings& settings) : beta_(Eigen::VectorXd::Zero(settings.get_K())), lambda_(-1), mean_(0) {}
  
  GAMParams(const Rcpp::List& state) { set_state(state); }
  Rcpp::List get_state() const {
    return Rcpp::List::create(_["beta"]=get_beta(), _["lambda"]=get_lambda(), _["mean"]=get_mean());
  }
  void set_state(const Rcpp::List& state) {
    set_beta(Rcpp::as<Eigen::VectorXd>(state["beta"]));
    set_lambda(Rcpp::as<double>(state["lambda"]));
    set_mean(Rcpp::as<double>(state["mean"]));
  }
  
  BINLESS_FORBID_COPY(GAMParams);
  
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, beta);
  BINLESS_GET_SET_DECL(double, double, lambda);
  BINLESS_GET_SET_DECL(double, double, mean);
  
};

// class that holds parameters that were output by fitting the mean of some statistic
struct MeanParams {
  template<typename Settings>
  MeanParams(const Settings& settings) : estimate_(Eigen::VectorXd::Zero(settings.get_nbins())), mean_(0) {}
  
  MeanParams(const Rcpp::List& state) { set_state(state); }
  Rcpp::List get_state() const {
    return Rcpp::List::create(_["estimate"]=get_estimate(), _["mean"]=get_mean());
  }
  void set_state(const Rcpp::List& state) {
    set_estimate(Rcpp::as<Eigen::VectorXd>(state["estimate"]));
    set_mean(Rcpp::as<double>(state["mean"]));
  }
  
  BINLESS_FORBID_COPY(MeanParams);
  
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, estimate);
  BINLESS_GET_SET_DECL(double, double, mean);
};


// Estimator is a policy class that takes data and some configuration info and performs a complete IRLS step on it
template<class SummarizerImpl>
class Summarizer : public SummarizerImpl {
public:
  template<typename FastData, typename Config>
  Summarizer(const FastData& data, const Config& conf) : SummarizerImpl(data,conf) {}
  
  //compute group sums of a vector of the size of the input data into the bins formed for the decay calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return SummarizerImpl::get_settings().get_binner()*vec; }
  
  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const std::vector<double>& log_expected, const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z, const Eigen::VectorXd& estimate);
  
};


template<class FitterImpl>
class Fitter : private FitterImpl {
public:
  template<typename FastData, typename Config>
  Fitter(const FastData& data, const Config& conf) : FitterImpl(data,conf) {}

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

template<class SummarizerImpl, class FitterImpl>
class Estimator : public Summarizer<SummarizerImpl>, public Fitter<FitterImpl> {
public:
  
  typedef Summarizer<SummarizerImpl> summarizer_t;
  typedef Fitter<FitterImpl> fitter_t;
  
  template<typename FastData, typename Config>
  Estimator(const FastData& data, const Config& conf) : summarizer_t(data,conf), fitter_t(data,conf) {}
  
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

#include "fast_estimator.ipp"

}
}

#endif

