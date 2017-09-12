#ifndef BOUNDS_CHECKER_HPP
#define BOUNDS_CHECKER_HPP

#include <Rcpp.h>

#include "BoundsOptimizer.hpp" //for bounds_t typedef
#include "Traits.hpp"
#include "BinnedData.hpp"
#include "util.hpp" // get_forbidden_values

//BoundsChecker does tag dispatching by type
template<typename> class BoundsChecker {};

//specialization in the case where the sign is positive
template<> class BoundsChecker<PositiveSign> {
public:
    BoundsChecker(const BinnedDataCore& binned) : minbeta_(Rcpp::min(binned.get_beta())) {};
    
    bool is_valid(bounds_t bounds) const {
        return bounds.first <= minbeta_;
    }

private:
    double minbeta_;
};

//specialization in the case where the sign is not constrained
template<> class BoundsChecker<AnySign> {
public:
    BoundsChecker(const BinnedDataCore&) {};
    bool is_valid(bounds_t) const { return true; }
};



//specialization in the case where degeneracies are forbidden
template<> class BoundsChecker<ForbidDegeneracy> {
public:
    BoundsChecker(const BinnedDataCore& binned) {
        Rcpp::NumericVector fv = get_forbidden_values(binned);
        minval_ = min(fv);
        maxval_ = max(fv);
        Rcpp::Rcout << " forbidden values: min= " << minval_ << " max= " << maxval_ << "\n";
    }
    
    bool is_valid(bounds_t bounds) const {
        return (bounds.first <= minval_) && (maxval_ <= bounds.second);
    }
    
private:
    double minval_, maxval_;
};

//specialization in the case where the degeneracies are allowed
template<> class BoundsChecker<AllowDegeneracy> {
public:
    BoundsChecker(const Rcpp::NumericVector&) {};
    bool is_valid(bounds_t) const { return true; }
};



#endif

