#ifndef FAST_DATAFRAME_HPP
#define FAST_DATAFRAME_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>

#include "Traits.hpp"
#include "FastData.hpp"
#include "fast_decay.hpp"
#include "fast_bias_mean.hpp"
#include "fast_exposure.hpp"

namespace binless {
namespace fast {

Rcpp::DataFrame get_as_dataframe(const FastData<Signal>& data, const ExposureEstimator& expo, const BiasEstimator& bias, const DecayEstimator& dec, const NumericVector lam1, double tol_val);

Rcpp::DataFrame get_as_dataframe(const FastData<Difference>& data, const NumericVector lam1, double tol_val);

}
}

#endif
