#ifndef FAST_BINLESS_HPP
#define FAST_BINLESS_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "fast_residuals.hpp"
#include "fast_distribution.hpp"

namespace binless {
namespace fast {
 
struct PrecisionPair { double abs,rel; };
struct SignalTriplet { std::vector<double> phihat, weights, beta; };
struct DifferenceQuadruplet { std::vector<double> deltahat,weights,delta,phi_ref; };

template<typename Lasso>
SignalTriplet step_signal(const FastSignalData& data, const ResidualsPair& z, std::vector<Lasso>& flo, const Eigen::ArrayXd& lam2, unsigned group=0);
template<typename Lasso>
DifferenceQuadruplet step_difference(const FastDifferenceData& data, const ResidualsPair& z, std::vector<Lasso>& flos, const Eigen::ArrayXd& lam2, unsigned ref);

PrecisionPair get_precision(const std::vector<double>& weights, const std::vector<double>& weights_old);
std::vector<double> remove_signal_degeneracy(const FastSignalData& data);
std::vector<double> shift_signal(const FastSignalData& data);

Rcpp::List binless(const DataFrame obs, unsigned nbins, const NumericVector lam2, double alpha, unsigned nouter=25, double tol_val=2e-1, unsigned bg_steps=5, unsigned free_decay=10000);

Rcpp::List binless_eval_cv(const List obs, const NumericVector lam2, double alpha, unsigned group=0, double tol_val=1e-1);

Rcpp::DataFrame binless_difference(const List obs, const NumericVector lam2, unsigned ref, double alpha, double tol_val=2e-1);

#include "fast_binless.ipp"

}
}

#endif

