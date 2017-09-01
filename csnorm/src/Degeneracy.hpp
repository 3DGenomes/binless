#ifndef DEGENERACY_HPP
#define DEGENERACY_HPP

#include <Rcpp.h>

class ForbidDegeneracy {
public:
    
    ForbidDegeneracy(const Rcpp::DataFrame& mat) : forbidden_values_(compute_forbidden_values(mat)) {}
    
    //exclude a number of patch values (used as UB borders) because they would make the system degenerate
    Rcpp::NumericVector filter_borders(const Rcpp::NumericVector& c) const;
    Rcpp::NumericVector get_forbidden_values() const { return forbidden_values_; }

private:
    //store a sorted list of values which have to be shrunk to zero
    Rcpp::NumericVector compute_forbidden_values(const Rcpp::DataFrame& mat) const;

    //return the minimum value encountered along each counter diagonal
    Rcpp::NumericVector get_minimum_diagonal_values(const Rcpp::NumericVector& value, const Rcpp::IntegerVector& diag_grp) const;

    Rcpp::NumericVector forbidden_values_;
};


class AllowDegeneracy {
public:
    Rcpp::NumericVector filter_borders(const Rcpp::NumericVector& c) const { return c; }
    Rcpp::NumericVector get_forbidden_values() { return Rcpp::NumericVector(); }
};


#endif

