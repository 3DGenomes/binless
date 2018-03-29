#ifndef FAST_ESTIMATOR_HPP
#define FAST_ESTIMATOR_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "spline.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {

class Estimator {
public:
  template<typename FastData, typename Config>
  Estimator(const FastData& data, const Config& conf) : 

};

}
}

#endif

