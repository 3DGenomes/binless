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

struct DecayEstimate {
  Eigen::VectorXd log_decay;
  DecaySummary summary;
  double lambda_diag;
};

DecayEstimate init_decay(const FastSignalData& data);

std::vector<double> compute_poisson_lsq_log_decay(const FastSignalData& data);

DecaySummary get_decay_summary(const FastSignalData& data, const DecayEstimate& dec);
void spline_log_decay_fit(const DecaySummary& summary, DecayEstimate& dec, double tol_val, const DecaySchedule& schedule);
void step_log_decay(const FastSignalData& data, DecayEstimate& dec, double tol_val);

}
}

#endif

