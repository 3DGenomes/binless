#ifndef FAST_DECAY_HPP
#define FAST_DECAY_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "spline.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {

struct DecayConfig {
  DecayConfig(double tol_val) : tol_val(tol_val) {}
  
  //parameters for decay calculation
  unsigned Kdiag=50;
  unsigned max_iter=100;
  double sigma=1;
  unsigned bins_per_bf=100;
  
  double tol_val; //default value will be overwritten
  
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
    //ncounts
    ncounts_ = binner_ * Eigen::VectorXd::Ones(log_distance_data.size());
    //compute mean log distance
    log_distance_ = ((binner_ * log_distance_data).array() / ncounts_.array()).matrix();
    //X: design matrix
    X_ = generate_spline_base(log_distance_, log_dmin_, log_dmax_, conf_.Kdiag);
    //D: build difference matrix
    D_ = second_order_difference_matrix(conf_.Kdiag);
    //C: build constraint matrix to forbid increase
    Cin_ = - first_order_difference_matrix(conf_.Kdiag);
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
  Eigen::VectorXd get_ncounts() const { return ncounts_; }
  Eigen::SparseMatrix<double> get_X() const { return X_; }
  Eigen::SparseMatrix<double> get_D() const { return D_; }
  Eigen::SparseMatrix<double> get_Cin() const { return Cin_; }
  
private:
  const DecayConfig& conf_;
  Eigen::VectorXd log_distance_;
  double log_dmin_, log_dmax_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd ncounts_;
  Eigen::SparseMatrix<double> X_,D_,Cin_; // design, difference and constraint matrices
};

struct DecaySummary {
  Eigen::VectorXd kappahat, weight;
};

struct DecayParams {
  DecayParams(const DecaySettings& settings) : beta_diag(Eigen::VectorXd::Zero(settings.get_Kdiag())),
    lambda_diag(-1), mean(0) {}
  Eigen::VectorXd beta_diag;
  double lambda_diag, mean;
};

class Decay {
public:
  
  //here, initialize flat log decay
  template<typename FastData>
  Decay(const FastData& data, const DecayConfig& conf) :
   settings_(data,conf), summary_(), params_(settings_),
   gam_(settings_.get_X(), settings_.get_D(), settings_.get_sigma())
    { gam_.set_inequality_constraints(settings_.get_Cin()); }
  
  
  //compute group sums of a vector of the size of the input data into the bins formed for the decay calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return settings_.get_binner()*vec; }
  
  Eigen::VectorXd get_beta_diag() const { return params_.beta_diag; }
  void set_beta_diag(const Eigen::VectorXd& beta) { params_.beta_diag = beta; }
  
  //get log decay along binned distances
  Eigen::VectorXd get_binned_log_decay() const {
    return settings_.get_X() * params_.beta_diag - Eigen::VectorXd::Constant(settings_.get_nbins(), params_.mean);
  }
  
  //get approximate log decay along distances in original data (same approx as during fitting)
  Eigen::VectorXd get_data_log_decay() const {
    return settings_.get_binner().transpose()*get_binned_log_decay();
  }
  
  double get_lambda_diag() const { return params_.lambda_diag; }
  void set_lambda_diag(double lambda_diag) { params_.lambda_diag = lambda_diag; }

  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const FastSignalData& data);
  //perform spline fit of summary data
  void update_params();
  //one complete IRLS iteration for log decay
  void step_log_decay(const FastSignalData& data) {
    update_summary(data);
    update_params();
  }
  
  
private:
  //compute average log decay (weighted by ncounts) in order to center it
  void center_log_decay() {
    params_.mean = settings_.get_ncounts().dot(settings_.get_X() * params_.beta_diag)/settings_.get_ncounts().sum();
  }
  
  const DecaySettings settings_; // parameters for performing the binning, constant
  DecaySummary summary_; // transformed data, iteration-specific
  DecayParams params_; // resulting fit, iteration-specific
  GeneralizedAdditiveModel gam_; //used to fit parameters
};

}
}

#endif

