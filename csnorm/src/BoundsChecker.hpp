#ifndef BOUNDS_CHECKER_HPP
#define BOUNDS_CHECKER_HPP

#include <Rcpp.h>

#include "BoundsComputer.hpp" //for bounds_t typedef
#include "Sign.hpp"

//BoundsChecker does tag dispatching by the type of Sign
template<typename Sign> class BoundsChecker {};

//specialization in the case where the sign is positive
template<> class BoundsChecker<PositiveSign> {
public:
    BoundsChecker(const Rcpp::NumericVector& beta) : beta_(beta) {};
    
    bool is_valid(bounds_t bounds) const {
        double minval = min(beta);
        return bounds.first <= minval;
    }

private:
    Rcpp::NumericVector beta_;
};

//specialization in the case where the sign is not constrained
template<> class BoundsChecker<AnySign> {
public:
    BoundsChecker(const Rcpp::NumericVector&) {};
    bool is_valid(bounds_t, const Rcpp::NumericVector&) const { return true; }
};


#endif

