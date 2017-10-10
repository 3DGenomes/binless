#ifndef FAST_BINLESS_HPP
#define FAST_BINLESS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"

struct ResidualsPair { std::vector<double> residuals,weights; };
struct PrecisionPair { double abs,rel; };
struct SignalTriplet { std::vector<double> phihat, weights, beta; };
struct DifferenceQuadruplet { std::vector<double> deltahat,weights,delta,phi_ref; };

template<typename FastData>
ResidualsPair get_normal_residuals(const FastData& data);
template<typename FastData>
ResidualsPair get_poisson_residuals(const FastData& data);

std::vector<double> fast_compute_poisson_exposures(const FastSignalData& data);
std::vector<double> fast_compute_exposures(const FastSignalData& data);
std::vector<double> fast_compute_log_biases(const FastSignalData& data);
std::vector<double> fast_compute_log_decay(const FastSignalData& data);
template<typename Lasso>
SignalTriplet fast_compute_signal(const FastSignalData& data, std::vector<Lasso>& flo, double lam2);
template<typename Lasso>
DifferenceQuadruplet fast_compute_difference(const FastDifferenceData& data, std::vector<Lasso>& flos, double lam2, unsigned ref);

PrecisionPair fast_precision(const std::vector<double>& weights, const std::vector<double>& weights_old);
std::vector<double> fast_remove_signal_degeneracy(const FastSignalData& data);
std::vector<double> fast_shift_signal(const FastSignalData& data);

Rcpp::List fast_binless(const DataFrame obs, unsigned nbins, unsigned nouter, double lam2, double tol_val);

Rcpp::DataFrame fast_binless_difference(const List obs, double lam2, double tol_val, unsigned ref);

#include "fast_binless.ipp"

#endif

