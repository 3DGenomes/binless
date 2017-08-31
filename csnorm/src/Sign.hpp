#ifndef SIGN_HPP
#define SIGN_HPP

#include <Rcpp.h>
#include <utility> //pair

struct PositiveSign {
    bool is_valid(std::pair<double,double> bounds, const Rcpp::NumericVector& beta) const {
        double minval = min(beta);
        return bounds.first <= minval;
    }
};

struct AnySign {
    bool is_valid(std::pair<double,double>, const Rcpp::NumericVector&) const { return true; }
};

#endif

