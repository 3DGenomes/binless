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


template<typename Leg>
class FitterImpl<Leg,Mean> {
public:
  template<typename SummarizerSettings, typename FastData, typename Config>
  FitterImpl(const SummarizerSettings& sset, const FastData& data, const Config& conf) : 
    settings_(sset, data, conf), params_(settings_) {}
  
  //update beta and lambda given phihat and weight
  void update_params(const Eigen::VectorXd& phihat, const Eigen::VectorXd& weight) {
    //center log_bias (no smoothing)
    double avg = get_settings().get_nobs().dot(phihat)/get_settings().get_nobs().sum();
    Eigen::VectorXd log_biases = (phihat.array() - avg).matrix();
    //cap estimates at 3SD from the mean
    double stdev = std::sqrt(log_biases.squaredNorm()/log_biases.rows());
    log_biases = log_biases.array().min(3*stdev).max(-3*stdev).matrix();
    
    set_estimate(log_biases);
  }
  
  //get X*beta
  Eigen::VectorXd get_estimate() const { return get_params().get_estimate(); }
  
  BINLESS_GET_CONSTREF_DECL(FitterSettings<Leg>, settings);
  BINLESS_GET_REF_DECL(Params<Mean>, params);
  
private:
  void set_estimate(const Eigen::VectorXd& estimate) { get_params().set_estimate(estimate); }
  
};

}
}

#endif

