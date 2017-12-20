#ifndef FAST_RESIDUALS_HPP
#define FAST_RESIDUALS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"

namespace binless {
namespace fast {

struct ResidualsPair { std::vector<double> residuals,weights; };

template<typename FastData>
ResidualsPair get_normal_residuals(const FastData& data);
template<typename FastData>
ResidualsPair get_poisson_residuals(const FastData& data);

#include "fast_residuals.ipp"

}
}

#endif

