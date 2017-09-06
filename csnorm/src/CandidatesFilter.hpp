#ifndef CANDIDATES_FILTER_HPP
#define CANDIDATES_FILTER_HPP

#include <Rcpp.h>

#include "BinnedData.hpp"
#include "Traits.hpp"

//class to filter upper bound candidates based on a degeneracy policy
//TODO: should the sign policy be included?
template<class Degeneracy> class CandidatesFilter {};

template<> class CandidatesFilter<AllowDegeneracy> {
public:
    CandidatesFilter(const BinnedData&) {}
    Rcpp::NumericVector filter(const Rcpp::NumericVector&) const { return Rcpp::NumericVector(); }
};


template<> class CandidatesFilter<ForbidDegeneracy> {
public:
    
    CandidatesFilter(const BinnedData& binned) : maxval_(max(get_forbidden_values(binned))) {}
    
    //exclude a number of patch values (used as UB borders) because they would make the system degenerate
    Rcpp::NumericVector filter(const Rcpp::NumericVector& candidates) const {
        std::vector<double> vc = Rcpp::as<std::vector<double> >(candidates);
        auto it = std::find_if(vc.begin(), vc.end(), [this](double d) { return d > maxval_; } );
        if (it == vc.begin() || it == vc.end()) Rcpp::stop("This should not happen because candidates should extend beyond patch values"); //TODO
        return Rcpp::wrap(std::vector<double>(it-1,vc.end()));
    }
    
private:
    double maxval_;
};


#endif

