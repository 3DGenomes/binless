#ifndef FAST_DECAY_HPP
#define FAST_DECAY_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "spline.hpp"

namespace binless {
namespace fast {

struct DecaySummary {
  Eigen::VectorXd distance, kappahat, weight, ncounts;
};

struct DecaySchedule {
  template<typename FastData>
  DecaySchedule(const FastData& data, double tol_val) {
    auto distance_std = data.get_distance();
    const Eigen::Map<const Eigen::VectorXd> distance(distance_std.data(),distance_std.size());
    Eigen::ArrayXd log_distance = distance.array().log();
    log_dmin = log_distance.minCoeff();
    log_dmax = log_distance.maxCoeff();
  }
  
  //parameters for decay calculation
  unsigned Kdiag=50;
  unsigned max_iter=100;
  double sigma=1;
  unsigned bins_per_bf=100;
  
  double log_dmin, log_dmax, tol_val; //default value will be overwritten
  
};

struct DecayParams {
  DecayParams(unsigned Kdiag) : beta_diag(Eigen::VectorXd::Zero(Kdiag)), lambda_diag(-1), mean(0) {}
  Eigen::VectorXd beta_diag;
  double lambda_diag, mean;
};

class Decay {
public:
  
  //here, initialize flat log decay
  template<typename FastData>
  Decay(const FastData& data, double tol_val) :
   schedule_(DecaySchedule(data, tol_val)), summary_(DecaySummary()), params_(DecayParams(schedule_.Kdiag)) {}
  
  Eigen::SparseMatrix<double> get_design(const Eigen::VectorXd& log_distance) const {
    return generate_spline_base(log_distance, schedule_.log_dmin, schedule_.log_dmax, schedule_.Kdiag);
  }
  
  Eigen::VectorXd get_beta_diag() const { return params_.beta_diag; }
  void set_beta_diag(const Eigen::VectorXd& beta) { params_.beta_diag = beta; }
  
  Eigen::VectorXd get_log_decay(const Eigen::VectorXd& log_distance) const {
    const Eigen::SparseMatrix<double> X = get_design(log_distance);
    auto log_decay = X * params_.beta_diag - Eigen::VectorXd::Constant(log_distance.rows(), params_.mean);
    return log_decay;
  }
  
  void center_log_decay(const Eigen::VectorXd& log_distance, const Eigen::VectorXd& w) {
    auto log_decay = get_log_decay(log_distance);
    double avg = w.dot(log_decay)/w.sum();
    //Rcpp::Rcout << "center log decay: before " << mean_ << " calc  " << avg << " after " << mean_+avg << "\n";
    params_.mean += avg;
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
  const DecaySchedule schedule_;
  DecaySummary summary_;
  DecayParams params_;
};

}
}

#endif

