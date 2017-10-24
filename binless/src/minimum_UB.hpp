#ifndef MINIMUM_UB_HPP
#define MINIMUM_UB_HPP

#include <Rcpp.h>
#include <vector>
#include <algorithm>

#include "CandidatesValues.hpp"
#include "Traits.hpp"
#include "BinnedData.hpp"

//get_minimum_UB_offset_sign
template<typename Offset, typename Sign> double get_minimum_UB_offset_sign(const BinnedDataCore& binned);

//get_minimum_UB_degen
double get_max_forbidden_border(const BinnedDataCore& binned, const Rcpp::NumericVector& beta);
template<typename Degeneracy> double get_minimum_UB_degen(const BinnedDataCore& binned, const Rcpp::NumericVector& borders);

//get_minimum_UB
template<typename Offset, //whether offset should be held at zero (ZeroOffset) or estimated (EstimatedOffset)
         typename Sign, //whether there is no constraint on the sign of the estimate (AnySign) or whether it must be positive (PositiveSign)
         typename Degeneracy> //whether to forbid certain values (ForbidDegeneracy) in order to avoid degeneracies in the model, or not (AllowDegeneracy)
double get_minimum_UB(const BinnedDataCore& binned) {
    const CandidatesValues<Offset,Sign> cv(binned);
    const Rcpp::NumericVector values = cv.get();
    const double minUBo = get_minimum_UB_offset_sign<Offset,Sign>(binned);
    const double minUBd = get_minimum_UB_degen<Degeneracy>(binned, values);
    const double minUB = std::max(minUBd,minUBo);
    //Rcpp::Rcout << "get_minimum UB: " << minUB << " = max(UBd= " << minUBd << ", UBo= " << minUBo << " )\n";
    return minUB;
}

#endif

