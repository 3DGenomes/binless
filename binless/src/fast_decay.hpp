#ifndef FAST_DECAY_HPP
#define FAST_DECAY_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"

namespace binless {
namespace fast {

struct DecaySummary {
  Eigen::VectorXd distance, kappahat, weight, ncounts;
};

struct DecaySchedule {
  unsigned Kdiag=50;
  unsigned max_iter=100;
  double sigma=1;
  //unsigned bins_per_bf=10;
};

class DecayEstimate {
public:
  DecayEstimate(const Eigen::VectorXd& log_decay, const DecaySummary& summary, double lambda_diag) :
    log_decay_(log_decay), summary_(summary), lambda_diag_(lambda_diag) {}
  
  Eigen::VectorXd get_log_decay(const Eigen::VectorXd&) const { return log_decay_; }
  void set_log_decay(const Eigen::VectorXd& log_decay) { log_decay_ = log_decay; }
  
  DecaySummary get_summary() const { return summary_; }
  void set_summary(const DecaySummary& summary) { summary_ = summary; }
  
  double get_lambda_diag() const { return lambda_diag_; }
  void set_lambda_diag(double lambda_diag) { lambda_diag_ = lambda_diag; }
  
private:
  Eigen::VectorXd log_decay_;
  DecaySummary summary_;
  double lambda_diag_;
};

DecayEstimate init_decay(unsigned nbins);

Eigen::VectorXd compute_poisson_lsq_log_decay(const FastSignalData& data, const DecayEstimate& dec);

DecaySummary get_decay_summary(const FastSignalData& data, const DecayEstimate& dec);
void spline_log_decay_fit(const DecaySummary& summary, DecayEstimate& dec, double tol_val, const DecaySchedule& schedule);
void step_log_decay(const FastSignalData& data, DecayEstimate& dec, double tol_val);

}
}

#endif

