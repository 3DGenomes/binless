#ifndef BOUNDS_CHECKER_HPP
#define BOUNDS_CHECKER_HPP

#include <Rcpp.h>

#include "BoundsComputer.hpp" //for bounds_t typedef
#include "Traits.hpp"
#include "Degeneracy.hpp"

//BoundsChecker does tag dispatching by type
template<typename> class BoundsChecker {};

//specialization in the case where the sign is positive
template<> class BoundsChecker<PositiveSign> {
public:
    BoundsChecker(const Rcpp::NumericVector& beta) : minbeta_(Rcpp::min(beta)) {};
    
    bool is_valid(bounds_t bounds) const {
        return bounds.first <= minbeta_;
    }

private:
    double minbeta_;
};

//specialization in the case where the sign is not constrained
template<> class BoundsChecker<AnySign> {
public:
    BoundsChecker(const Rcpp::NumericVector&) {};
    bool is_valid(bounds_t) const { return true; }
};



//specialization in the case where degeneracies are forbidden
template<> class BoundsChecker<ForbidDegeneracy> {
public:
    BoundsChecker(const Rcpp::NumericVector& forbidden_values) :
     minval_(min(forbidden_values)), maxval_(max(forbidden_values)) {};
    
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

