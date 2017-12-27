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
  unsigned Kdiag=50;
  unsigned max_iter=100;
  double sigma=1;
  unsigned bins_per_bf=100;
};

class DecayEstimate {
public:
  DecayEstimate(double log_dmin, double log_dmax, unsigned Kdiag, const DecaySummary& summary, double lambda_diag) :
    log_dmin_(log_dmin), log_dmax_(log_dmax), Kdiag_(Kdiag), beta_diag_(Eigen::VectorXd::Zero(Kdiag_)),
    summary_(summary), lambda_diag_(lambda_diag), mean_(0) {}
  
  Eigen::SparseMatrix<double> get_design(const Eigen::VectorXd& log_distance) const {
    return generate_spline_base(log_distance, log_dmin_, log_dmax_, Kdiag_);
  }
  
  Eigen::VectorXd get_beta_diag() const { return beta_diag_; }
  void set_beta_diag(const Eigen::VectorXd& beta) {
    beta_diag_ = beta;
  }
  
  Eigen::VectorXd get_log_decay(const Eigen::VectorXd& log_distance) const {
    const Eigen::SparseMatrix<double> X = get_design(log_distance);
    auto log_decay = X * beta_diag_ - Eigen::VectorXd::Constant(log_distance.rows(),mean_);
    return log_decay;
  }
  
  void center_log_decay(const Eigen::VectorXd& log_distance, const Eigen::VectorXd& w) {
    auto log_decay = get_log_decay(log_distance);
    double avg = w.dot(log_decay)/w.sum();
    //Rcpp::Rcout << "center log decay: before " << mean_ << " calc  " << avg << " after " << mean_+avg << "\n";
    mean_ += avg;
  }
  
  DecaySummary get_summary() const { return summary_; }
  void set_summary(const DecaySummary& summary) { summary_ = summary; }
  
  double get_lambda_diag() const { return lambda_diag_; }
  void set_lambda_diag(double lambda_diag) { lambda_diag_ = lambda_diag; }
  
private:
  double log_dmin_, log_dmax_;
  unsigned Kdiag_;
  Eigen::VectorXd beta_diag_;
  DecaySummary summary_;
  double lambda_diag_, mean_;
};

//here, initialize flat log decay
template<typename FastData>
DecayEstimate init_decay(const FastData& data) {
  DecaySchedule schedule;
  DecaySummary summary;
  auto distance_std = data.get_distance();
  const Eigen::Map<const Eigen::VectorXd> distance(distance_std.data(),distance_std.size());
  Eigen::ArrayXd log_distance = distance.array().log();
  double ldmin = log_distance.minCoeff();
  double ldmax = log_distance.maxCoeff();
  DecayEstimate dec(ldmin, ldmax, schedule.Kdiag, summary, -1);
  return dec;
}

DecaySummary compute_poisson_lsq_log_decay(const FastSignalData& data, const DecayEstimate& dec, const DecaySchedule& schedule);

DecaySummary get_decay_summary(const FastSignalData& data, const DecayEstimate& dec, const DecaySchedule& schedule);
void spline_log_decay_fit(const DecaySummary& summary, DecayEstimate& dec, double tol_val, const DecaySchedule& schedule);
void step_log_decay(const FastSignalData& data, DecayEstimate& dec, double tol_val);

}
}

#endif

