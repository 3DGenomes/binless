#ifndef FAST_BINLESS_HPP
#define FAST_BINLESS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"

struct ResidualsPair { std::vector<double> residuals,weights; };
struct SignalTriplet { std::vector<double> phihat, weights, beta; };

ResidualsPair get_normal_residuals(const FastData& data);
ResidualsPair get_poisson_residuals(const FastData& data);
std::vector<double> fast_compute_exposures(const FastData& data);
std::vector<double> fast_compute_log_biases(const FastData& data);
std::vector<double> fast_compute_log_decay(const FastData& data);
template<typename Lasso>
SignalTriplet fast_compute_signal(const FastData& data, std::vector<Lasso>& flo, double lam2);
double fast_precision(const std::vector<double>& weights, const std::vector<double>& weights_old);
std::vector<double> fast_remove_signal_degeneracy(const FastData& data);
std::vector<double> fast_shift_signal(const FastData& data);

List fast_binless(const DataFrame obs, unsigned nbins, unsigned nouter, double lam2, double tol_val);

#include "fast_binless.ipp"

#endif

