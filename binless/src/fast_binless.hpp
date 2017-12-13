#ifndef FAST_BINLESS_HPP
#define FAST_BINLESS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"

namespace binless {
namespace fast {

struct ResidualsPair { std::vector<double> residuals,weights; };
struct PrecisionPair { double abs,rel; };
struct SignalTriplet { std::vector<double> phihat, weights, beta; };
struct DifferenceQuadruplet { std::vector<double> deltahat,weights,delta,phi_ref; };
struct DecaySummary { std::vector<double> distance, kappahat, weight; };
struct DecayFit {
  std::vector<double> log_decay;
  DecaySummary dec;
  double lambda_diag;
};

template<typename FastData>
ResidualsPair get_normal_residuals(const FastData& data);
template<typename FastData>
ResidualsPair get_poisson_residuals(const FastData& data);

std::vector<double> compute_poisson_lsq_exposures(const FastSignalData& data);
std::vector<double> step_exposures(const FastSignalData& data);
std::vector<double> compute_poisson_lsq_log_biases(const FastSignalData& data);
std::vector<double> step_log_biases(const FastSignalData& data);
std::vector<double> compute_poisson_lsq_log_decay(const FastSignalData& data);
DecaySummary get_decay_summary(const FastSignalData& data);
DecayFit pointwise_log_decay_fit(const DecaySummary& dec);
std::vector<double> step_log_decay(const FastSignalData& data);
template<typename Lasso>
SignalTriplet step_signal(const FastSignalData& data, std::vector<Lasso>& flo, double lam2, unsigned group=0);
template<typename Lasso>
DifferenceQuadruplet step_difference(const FastDifferenceData& data, std::vector<Lasso>& flos, double lam2, unsigned ref);

PrecisionPair get_precision(const std::vector<double>& weights, const std::vector<double>& weights_old);
std::vector<double> remove_signal_degeneracy(const FastSignalData& data);
std::vector<double> shift_signal(const FastSignalData& data);

Rcpp::List binless(const DataFrame obs, unsigned nbins, double lam2, unsigned nouter=20, double tol_val=1e-1, unsigned bg_steps=5);

Rcpp::List binless_eval_cv(const List obs, const NumericVector lam2, unsigned group=0, double tol_val=1e-1);

Rcpp::DataFrame binless_difference(const List obs, double lam2, unsigned ref, double tol_val=1e-1);

#include "fast_binless.ipp"

}
}

#endif

