#ifndef FAST_FITTER_MEAN_HPP
#define FAST_FITTER_MEAN_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "macros.hpp"
#include "Traits.hpp"

namespace binless {
namespace fast {

// class that holds parameters that were output by fitting the mean of some statistic
template<>
struct Params<Mean> {
  template<typename Settings>
  Params(const Settings& settings) : estimate_(Eigen::VectorXd::Zero(settings.get_nbins())), mean_(0) {}
  
  Params(const Rcpp::List& state) { set_state(state); }
  Rcpp::List get_state() const {
    return Rcpp::List::create(_["estimate"]=get_estimate(), _["mean"]=get_mean());
  }
  void set_state(const Rcpp::List& state) {
    set_estimate(Rcpp::as<Eigen::VectorXd>(state["estimate"]));
    set_mean(Rcpp::as<double>(state["mean"]));
  }
  
  BINLESS_FORBID_COPY(Params);
  
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, estimate);
  BINLESS_GET_SET_DECL(double, double, mean);
};

// class that holds the data used to fit summaries and generate corresponding params
// cannot be built directly, as it is meant to be constructed by each
// domain-specific child
template<>
class FitterSettings<Mean> {
  
  BINLESS_GET_CONSTREF_DECL(unsigned, nbins);
  BINLESS_GET_CONSTREF_DECL(Eigen::VectorXd, nobs);
  
protected:
  FitterSettings(unsigned nbins, const Eigen::VectorXd& nobs) : nbins_(nbins), nobs_(nobs) {}
};

template<typename Leg>
class FitterImpl<Leg,Mean> {
public:
  FitterImpl(const SummarizerSettings& sset, const Config<Leg,Mean>& conf) : 
    settings_(FitterSettingsImpl<Leg,Mean>(sset, conf)), params_(settings_) {}
  
  //update beta and lambda given phihat and weight
  void update_params(const Eigen::VectorXd& phihat, const Eigen::VectorXd& weight) {
    //the mean is actually already calculated: phihat
    Eigen::VectorXd estimate = phihat;
    
    //center estimate based on number of observations in data
    if (FitterTraits<Leg,Mean>::center) {
      double avg = get_settings().get_nobs().dot(estimate)/get_settings().get_nobs().sum();
      estimate = (estimate.array() - avg).matrix();
    }
    //cap estimates at 3SD from the mean
    if (FitterTraits<Leg,Mean>::cap) {
      double stdev = std::sqrt(estimate.squaredNorm()/estimate.rows());
      estimate = estimate.array().min(3*stdev).max(-3*stdev).matrix();
    }
    
    //print debug info
    if (FitterTraits<Leg,Mean>::debug) {
      Rcpp::Rcout << "phihat: " << phihat.transpose() << "\n";
      Rcpp::Rcout << "weight: " << weight.transpose() << "\n";
      Rcpp::Rcout << "estimate: " << estimate.transpose() << "\n";
    }
    //this estimate is never centered after fitting (internal centering)
    set_estimate(estimate);
  }
  
  //get X*beta
  Eigen::VectorXd get_estimate() const { return get_params().get_estimate(); }
  
  BINLESS_GET_CONSTREF_DECL(FitterSettings<Mean>, settings);
  BINLESS_GET_REF_DECL(Params<Mean>, params);
  
private:
  void set_estimate(const Eigen::VectorXd& estimate) { get_params().set_estimate(estimate); }
  
};

}
}

#endif

