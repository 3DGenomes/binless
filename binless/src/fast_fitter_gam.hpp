#ifndef FAST_FITTER_GAM_HPP
#define FAST_FITTER_GAM_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "macros.hpp"
#include "gam.hpp"
#include "Traits.hpp"

namespace binless {
namespace fast {


// class that holds parameters that were output by a GAM fit
template<>
struct Params<GAM> {
  template<typename Settings>
  Params(const Settings& settings) : beta_(Eigen::VectorXd::Zero(settings.get_K())), lambda_(-1), mean_(0) {}
  
  Params(const Rcpp::List& state) { set_state(state); }
  Rcpp::List get_state() const {
    return Rcpp::List::create(_["beta"]=get_beta(), _["lambda"]=get_lambda(), _["mean"]=get_mean());
  }
  void set_state(const Rcpp::List& state) {
    set_beta(Rcpp::as<Eigen::VectorXd>(state["beta"]));
    set_lambda(Rcpp::as<double>(state["lambda"]));
    set_mean(Rcpp::as<double>(state["mean"]));
  }
  
  BINLESS_FORBID_COPY(Params);
  
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, beta);
  BINLESS_GET_SET_DECL(double, double, lambda);
  BINLESS_GET_SET_DECL(double, double, mean);
  
};

// class that holds the data used to fit summaries and generate corresponding params
// cannot be built directly, as it is meant to be constructed by each
// domain-specific child
template<>
class FitterSettings<GAM> {
  
  BINLESS_GET_CONSTREF_DECL(unsigned, max_iter);
  BINLESS_GET_CONSTREF_DECL(double, tol_val);
  BINLESS_GET_CONSTREF_DECL(double, sigma);
  BINLESS_GET_CONSTREF_DECL(double, K);
  BINLESS_GET_CONSTREF_DECL(unsigned, nbins);
  BINLESS_GET_CONSTREF_DECL(Eigen::VectorXd, nobs);
  
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, X);
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, D);
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, Cin);
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, Ceq);
  
  
protected:
  FitterSettings(unsigned max_iter, double tol_val, double sigma, unsigned K, unsigned nbins, const Eigen::VectorXd& nobs) :
    max_iter_(max_iter), tol_val_(tol_val), sigma_(sigma), K_(K), nbins_(nbins), nobs_(nobs) {}
};

template<typename Leg>
class FitterImpl<Leg,GAM> {
public:
  FitterImpl(const SummarizerSettings& sset, const Config<Leg,GAM>& conf) : 
    settings_(FitterSettingsImpl<Leg,GAM>(sset, conf)), params_(settings_), gam_(settings_.get_X(), settings_.get_D(), settings_.get_sigma())
  {
    if (FitterTraits<Leg,GAM>::has_inequality_constraints) gam_.set_inequality_constraints(get_settings().get_Cin());
    if (FitterTraits<Leg,GAM>::has_equality_constraints) gam_.set_equality_constraints(get_settings().get_Ceq());
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
    //this estimate is always centered
    center_estimate(); 
  }
  
  //get X*beta
  Eigen::VectorXd get_estimate() const { return get_settings().get_X() * get_params().get_beta(); }
  
  BINLESS_GET_CONSTREF_DECL(FitterSettings<GAM>, settings);
  BINLESS_GET_REF_DECL(Params<GAM>, params);
  
private:
  Eigen::VectorXd get_beta() const { return get_params().get_beta(); }
  void set_beta(const Eigen::VectorXd& beta) { get_params().set_beta(beta); }
  
  double get_lambda() const { return get_params().get_lambda(); }
  void set_lambda(double lambda) { get_params().set_lambda(lambda); }
  
  //compute average of estimate in order to center it
  void center_estimate() {
    get_params().set_mean(get_settings().get_nobs().dot(get_estimate())
                          /get_settings().get_nobs().sum());
  }
  
private:
  GeneralizedAdditiveModel<typename FitterTraits<Leg,GAM>::library> gam_; //used to fit parameters
};

}
}

#endif

