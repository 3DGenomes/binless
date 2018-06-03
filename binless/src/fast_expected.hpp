#ifndef FAST_EXPECTED_HPP
#define FAST_EXPECTED_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Traits.hpp"
#include "fast_decay.hpp"
#include "fast_bias_mean.hpp"
#include "fast_exposure.hpp"

namespace binless {
namespace fast {

// returns the expected log mean for each observation (must be multiplied by nobs to match with observed counts per bin)
template<typename Derived>
std::vector<double> get_log_expected(const FastData<Derived>& data, const ExposureEstimator& expo, const BiasEstimator& bias, const DecayEstimator& dec);

#include "fast_expected.ipp"

}
}

#endif

