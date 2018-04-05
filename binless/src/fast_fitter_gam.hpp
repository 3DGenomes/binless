#ifndef FAST_FITTER_GAM_HPP
#define FAST_FITTER_GAM_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "macros.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {

template<typename Leg>
struct GAMFitterTraits;

template<typename Leg>
class FitterSettings;

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

template<typename Leg>
class GAMFitterImpl {
public:
  template<typename SummarizerSettings, typename FastData, typename Config>
  GAMFitterImpl(const SummarizerSettings& sset, const FastData& data, const Config& conf) : 
    settings_(sset, data, conf), params_(settings_), gam_(settings_.get_X(), settings_.get_D(), settings_.get_sigma())
  {
    if (GAMFitterTraits<Leg>::has_inequality_constraints) gam_.set_inequality_constraints(get_settings().get_Cin());
    if (GAMFitterTraits<Leg>::has_equality_constraints) gam_.set_equality_constraints(get_settings().get_Ceq());
  }
  
  //update beta and lambda given phihat and weight
  void update_params(const Eigen::VectorXd& phihat, const Eigen::VectorXd& weight) {
    //extract data
    const Eigen::VectorXd y(phihat);
    const Eigen::VectorXd Sm1 = weight.array().sqrt().matrix();
    gam_.optimize(y,Sm1,settings_.get_max_iter(), settings_.get_tol_val());
    //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
    set_lambda(gam_.get_lambda());
    set_beta(gam_.get_beta());
  }
  
  //get X*beta
  Eigen::VectorXd get_estimate() const { return get_settings().get_X() * get_params().get_beta(); }
  
  BINLESS_GET_CONSTREF_DECL(FitterSettings<Leg>, settings);
  BINLESS_GET_REF_DECL(GAMParams, params);
  
private:
  Eigen::VectorXd get_beta() const { return get_params().get_beta(); }
  void set_beta(const Eigen::VectorXd& beta) { get_params().set_beta(beta); }
  
  double get_lambda() const { return get_params().get_lambda(); }
  void set_lambda(double lambda) { get_params().set_lambda(lambda); }
  
private:
  GeneralizedAdditiveModel<typename GAMFitterTraits<Leg>::library> gam_; //used to fit parameters
};

}
}

#endif

