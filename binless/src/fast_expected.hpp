#ifndef FAST_EXPECTED_HPP
#define FAST_EXPECTED_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Traits.hpp"
#include "fast_decay.hpp"

namespace binless {
namespace fast {

template<typename Derived>
std::vector<double> get_log_expected(const FastData<Derived>& data, const DecayEstimator& dec);

#include "fast_expected.ipp"

}
}

#endif

