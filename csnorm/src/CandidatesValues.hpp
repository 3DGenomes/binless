#ifndef CANDIDATES_VALUES_HPP
#define CANDIDATES_VALUES_HPP

#include <Rcpp.h>

#include "Traits.hpp"
#include "BinnedData.hpp"

template<typename Offset, typename Sign> class CandidatesValues {};

//when offset=0, UB borders are on the folded scale
template<typename Sign> class CandidatesValues<ZeroOffset, Sign> {
public:
    CandidatesValues(const BinnedDataCore& binned) : binned_(binned) {}
    Rcpp::NumericVector get() const { return abs(binned_.get_beta()); }
private:
    const BinnedDataCore& binned_;
};


//when sign>0 and the offset is estimated, UB borders are on the natural scale
template<> class CandidatesValues<EstimatedOffset, PositiveSign> {
public:
    CandidatesValues(const BinnedDataCore& binned) : binned_(binned) {}
    Rcpp::NumericVector get() const { return binned_.get_beta(); }
private:
    const BinnedDataCore& binned_;
};



#endif

