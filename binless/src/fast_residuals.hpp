#ifndef FAST_RESIDUALS_HPP
#define FAST_RESIDUALS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "fast_expected.hpp"

namespace binless {
namespace fast {

struct ResidualsPair { std::vector<double> residuals,weights; };

template<typename FastData>
ResidualsPair get_normal_residuals(const FastData& data, const DecayEstimator& dec);
template<typename FastData>
ResidualsPair get_poisson_residuals(const FastData& data, const DecayEstimator& dec);

#include "fast_residuals.ipp"

}
}

#endif

