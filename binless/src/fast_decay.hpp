#ifndef FAST_DECAY_HPP
#define FAST_DECAY_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"

namespace binless {
namespace fast {

struct DecaySummary { std::vector<double> distance, kappahat, weight; };
struct DecayFit {
  std::vector<double> log_decay;
  DecaySummary dec;
  double lambda_diag;
};

std::vector<double> compute_poisson_lsq_log_decay(const FastSignalData& data);

DecaySummary get_decay_summary(const FastSignalData& data);
DecayFit spline_log_decay_fit(const DecaySummary& dec, double tol_val, unsigned Kdiag=50, unsigned max_iter=100, double sigma=0.1);
std::vector<double> step_log_decay(const FastSignalData& data, double tol_val);

}
}

#endif

