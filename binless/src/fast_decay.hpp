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

class DecaySchedule {
  
public:
  template<typename FastData>
  DecaySchedule(const FastData& data, const DecayConfig& conf) : conf_(conf) {
    //log_distance and bounds
    auto distance_std = data.get_distance();
    const Eigen::Map<const Eigen::VectorXd> distance_data(distance_std.data(),distance_std.size());
    Eigen::VectorXd log_distance_data = distance_data.array().log();
    log_dmin_ = log_distance_data.minCoeff();
    log_dmax_ = log_distance_data.maxCoeff();
    //binner matrix
    binner_ = bin_data_evenly(log_distance_data, get_nbins(), true); //true to drop unused bins
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
  unsigned get_nbins() const { return conf_.Kdiag*conf_.bins_per_bf; }
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
  Eigen::VectorXd ncounts_;
  Eigen::SparseMatrix<double> X_,D_,Cin_; // design, difference and constraint matrices
};

struct DecaySummary {
  Eigen::VectorXd kappahat, weight;
};

struct DecayParams {
  DecayParams(const DecaySchedule& schedule) : beta_diag(Eigen::VectorXd::Zero(schedule.get_Kdiag())),
    lambda_diag(-1), mean(0) {}
  Eigen::VectorXd beta_diag;
  double lambda_diag, mean;
};

class Decay {
public:
  
  //here, initialize flat log decay
  template<typename FastData>
  Decay(const FastData& data, const DecayConfig& conf) :
   schedule_(DecaySchedule(data,conf)), summary_(DecaySummary()), params_(DecayParams(schedule_)) {}
  
  //compute group sums of a vector of the size of the input data into the bins formed for the decay calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return schedule_.get_binner()*vec; }
  
  Eigen::SparseMatrix<double> get_design(const Eigen::VectorXd& log_distance) const {
    return generate_spline_base(log_distance, schedule_.get_log_dmin(), schedule_.get_log_dmax(), schedule_.get_Kdiag());
  }
  
  Eigen::VectorXd get_beta_diag() const { return params_.beta_diag; }
  void set_beta_diag(const Eigen::VectorXd& beta) { params_.beta_diag = beta; }
  
  Eigen::VectorXd get_log_decay(const Eigen::VectorXd& log_distance) const {
    const Eigen::SparseMatrix<double> X = get_design(log_distance);
    auto log_decay = X * params_.beta_diag - Eigen::VectorXd::Constant(log_distance.rows(), params_.mean);
    return log_decay;
  }
  
  DecaySummary get_summary() const { return summary_; }
  void set_summary(const DecaySummary& summary) { summary_ = summary; }
  
  double get_lambda_diag() const { return params_.lambda_diag; }
  void set_lambda_diag(double lambda_diag) { params_.lambda_diag = lambda_diag; }

  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const FastSignalData& data);
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
  const DecaySchedule schedule_; // parameters for performing the binning, constant, data-independent
  DecaySummary summary_; // transformed data, iteration-specific
  DecayParams params_; // resulting fit, iteration-specific
};

}
}

#endif

