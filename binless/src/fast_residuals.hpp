#ifndef FAST_RESIDUALS_HPP
#define FAST_RESIDUALS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "fast_expected.hpp"
#include "fast_distribution.hpp"

namespace binless {
namespace fast {

struct ResidualsPair { std::vector<double> residuals,weights; };

//residuals: normal with log-link, 0 drops data
template<typename FastData>
ResidualsPair get_residuals(const NormalDistribution& dist, const FastData& data, const DecayEstimator& dec);

//residuals: poisson with log-link
template<typename FastData>
ResidualsPair get_residuals(const PoissonDistribution& dist, const FastData& data, const DecayEstimator& dec);

//residuals: negative binomial  with log-link
template<typename FastData>
ResidualsPair get_residuals(const NegativeBinomialDistribution& dist, const FastData& data, const DecayEstimator& dec);

#include "fast_residuals.ipp"

}
}

#endif

