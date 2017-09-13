#ifndef CANDIDATES_UB_MIN_OSI_HPP
#define CANDIDATES_UB_MIN_OSI_HPP

#include <Rcpp.h>
#include <algorithm>

#include "Traits.hpp"
#include "BinnedData.hpp"

template<typename Offset, typename Sign> class CandidatesUBMinOsi {};

//offset=0 and sign>0
template<> class CandidatesUBMinOsi<ZeroOffset, PositiveSign> {
public:
    CandidatesUBMinOsi(const BinnedDataCore& binned) : binned_(binned) {}
    double get() const { return std::max( - (double)(Rcpp::min(binned_.get_beta())), 0.); }
private:
    const BinnedDataCore& binned_;
};

//offset=0 and any sign
template<> class CandidatesUBMinOsi<ZeroOffset, AnySign> {
public:
    CandidatesUBMinOsi(const BinnedDataCore&) {}
    double get() const { return 0; }
};

//offset estimated and sign>0
template<> class CandidatesUBMinOsi<EstimatedOffset, PositiveSign> {
public:
    CandidatesUBMinOsi(const BinnedDataCore& binned) : binned_(binned) {}
    double get() const { return Rcpp::min(binned_.get_beta()); }
private:
    const BinnedDataCore& binned_;
};


#endif

