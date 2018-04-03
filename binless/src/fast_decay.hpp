#ifndef FAST_DECAY_HPP
#define FAST_DECAY_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "spline.hpp"
#include "gam.hpp"
#include "QuadProgGAMLibrary.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

struct DecayConfig {
  DecayConfig(double tol_val, double free_decay) : tol_val(tol_val), free_decay(free_decay) {}
  
  //parameters for decay calculation
  unsigned Kdiag=50;
  unsigned max_iter=100;
  double sigma=1;
  unsigned bins_per_bf=100;
  
  //default values will be overwritten
  double tol_val;      // tolerance on decay convergence
  unsigned free_decay; // distance in bases until which decay is not forced to decrease
  
};

class DecaySettings {
  
public:
  template<typename FastData>
  DecaySettings(const FastData& data, const DecayConfig& conf) : conf_(conf) {
    //log_distance and bounds
    auto distance_std = data.get_distance();
    const Eigen::Map<const Eigen::VectorXd> distance_data(distance_std.data(),distance_std.size());
    Eigen::VectorXd log_distance_data = distance_data.array().log();
    log_dmin_ = log_distance_data.minCoeff();
    log_dmax_ = log_distance_data.maxCoeff();
    //binner matrix
    binner_ = bin_data_evenly(log_distance_data, conf_.Kdiag*conf_.bins_per_bf, true); //true to drop unused bins
    nbins_ = binner_.rows();
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    nobs_ = binner_ * nobs_data;
    //compute mean log distance
    log_distance_ = ((binner_ * (log_distance_data.array() * nobs_data.array()).matrix()).array() / nobs_.array()).matrix();
    //X: design matrix
    X_ = generate_spline_base(log_distance_, log_dmin_, log_dmax_, conf_.Kdiag);
    //D: build difference matrix
    D_ = second_order_difference_matrix(conf_.Kdiag);
    //C: build constraint matrix to forbid increase
    unsigned free_first = conf_.Kdiag * (std::log(conf_.free_decay)-log_dmin_)/(log_dmax_-log_dmin_);
    //Rcpp::Rcout << "Free decay in " << free_first << " out of " << conf_.Kdiag << " basis functions\n";
    Cin_ = decreasing_constraint(conf_.Kdiag, free_first);
  }
  
  double get_Kdiag() const { return conf_.Kdiag; }
  unsigned get_nbins() const { return nbins_; }
  unsigned get_max_iter() const { return conf_.max_iter; }
  double get_tol_val() const { return conf_.tol_val; }
  double get_sigma() const { return conf_.sigma; }
  
  Eigen::VectorXd get_log_distance() const { return log_distance_; }
  double get_log_dmin() const { return log_dmin_; }
  double get_log_dmax() const { return log_dmax_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  Eigen::SparseMatrix<double> get_X() const { return X_; }
  Eigen::SparseMatrix<double> get_D() const { return D_; }
  Eigen::SparseMatrix<double> get_Cin() const { return Cin_; }
  
private:
  const DecayConfig& conf_;
  Eigen::VectorXd log_distance_;
  double log_dmin_, log_dmax_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd nobs_;
  Eigen::SparseMatrix<double> X_,D_,Cin_; // design, difference and constraint matrices
};

struct DecaySummary {
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, kappahat);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, weight);
};

struct DecayParams {
  
  DecayParams(const DecaySettings& settings) : beta_(Eigen::VectorXd::Zero(settings.get_Kdiag())), lambda_(-1), mean_(0) {}

  DecayParams(const Rcpp::List& state) { set_state(state); }
  Rcpp::List get_state() const {
    return Rcpp::List::create(_["beta"]=get_beta(), _["lambda"]=get_lambda(), _["mean"]=get_mean());
  }
  void set_state(const Rcpp::List& state) {
    set_beta(Rcpp::as<Eigen::VectorXd>(state["beta"]));
    set_lambda(Rcpp::as<double>(state["lambda"]));
    set_mean(Rcpp::as<double>(state["mean"]));
  }

  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, beta);
  BINLESS_GET_SET_DECL(double, double, lambda);
  BINLESS_GET_SET_DECL(double, double, mean);
  
};

class DecayEstimator {
public:
  
  //here, initialize flat log decay
  template<typename FastData>
  DecayEstimator(const FastData& data, const DecayConfig& conf) :
   settings_(data,conf), summary_(), params_(settings_),
   gam_(settings_.get_X(), settings_.get_D(), settings_.get_sigma())
    { gam_.set_inequality_constraints(settings_.get_Cin()); }
  
  
  //compute group sums of a vector of the size of the input data into the bins formed for the decay calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return settings_.get_binner()*vec; }
   
  //get log decay along binned distances
  Eigen::VectorXd get_binned_estimate() const {
    return settings_.get_X() * params_.get_beta() - Eigen::VectorXd::Constant(settings_.get_nbins(), params_.get_mean());
  }
  
  //get approximate log decay along distances in original data (same approx as during fitting)
  Eigen::VectorXd get_data_estimate() const {
    return settings_.get_binner().transpose()*get_binned_estimate();
  }
  
  //provide a way to store and recall the state of the estimator
  Rcpp::List get_state() const { return params_.get_state(); }
  void set_state(const Rcpp::List& state) { params_.set_state(state); }

  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const std::vector<double>& log_expected, const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z);
  //perform spline fit of summary data
  void update_params();
  
private:
  Eigen::VectorXd get_beta() const { return params_.get_beta(); }
  void set_beta(const Eigen::VectorXd& beta) { params_.set_beta(beta); }
  
  double get_lambda() const { return params_.get_lambda(); }
  void set_lambda(double lambda) { params_.set_lambda(lambda); }

  //compute average log decay (weighted by nobs) in order to center it
  void center_estimate() {
    params_.set_mean(settings_.get_nobs().dot(settings_.get_X() * params_.get_beta())/settings_.get_nobs().sum());
  }
  
  const DecaySettings settings_; // parameters for performing the binning, constant
  DecaySummary summary_; // transformed data, iteration-specific
  DecayParams params_; // resulting fit, iteration-specific
  GeneralizedAdditiveModel<QuadProgGAMLibrary> gam_; //used to fit parameters
};

}
}

#endif

