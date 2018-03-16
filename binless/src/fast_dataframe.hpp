#ifndef FAST_DATAFRAME_HPP
#define FAST_DATAFRAME_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>

#include "Traits.hpp"
#include "FastData.hpp"
#include "fast_decay.hpp"

namespace binless {
namespace fast {

Rcpp::DataFrame get_as_dataframe(const FastData<Signal>& data, const DecayEstimator& dec, double tol_val);

Rcpp::DataFrame get_as_dataframe(const FastData<Difference>& data, double tol_val);

}
}

#endif
