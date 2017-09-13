#ifndef CANDIDATES_UB_MIN_DEGEN_HPP
#define CANDIDATES_UB_MIN_DEGEN_HPP

#include <Rcpp.h>

#include "Traits.hpp"
#include "BinnedData.hpp"

double get_max_forbidden_border(const BinnedDataCore& binned, const Rcpp::NumericVector& beta);

template<typename Degeneracy> class CandidatesUBMinDegen {};

template<> class CandidatesUBMinDegen<AllowDegeneracy> {
public:
    CandidatesUBMinDegen(const BinnedDataCore&) {}
    double get(const Rcpp::NumericVector& borders) const { return min(borders)-1; }
};


template<> class CandidatesUBMinDegen<ForbidDegeneracy> {
public:
    CandidatesUBMinDegen(const BinnedDataCore& binned) : binned_(binned) {}
    double get(const Rcpp::NumericVector& borders) const {
        return get_max_forbidden_border(binned_, borders); //folded or not, depending on borders
    }
private:
    
    const BinnedDataCore& binned_;
};

#endif

