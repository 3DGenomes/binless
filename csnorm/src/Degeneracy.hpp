#ifndef DEGENERACY_HPP
#define DEGENERACY_HPP

#include <Rcpp.h>
#include "util.hpp"

struct ForbidDegeneracy {
    //return a sorted list of values which have to be shrunk to zero
    Rcpp::NumericVector get_forbidden_values(const DataFrame& mat) {
        Rcpp::NumericVector beta = mat["beta"];
        Rcpp::IntegerVector diag_grp = mat["diag.grp"];
        return get_minimum_diagonal_values(beta, diag_grp); //TODO: beta or abs(beta) or a combination?
    }
    
    //exclude a number of patch values (used as UB borders) because they would make the system degenerate
    Rcpp::NumericVector filter_borders(const Rcpp::NumericVector& c, const Rcpp::NumericVector& fv) {
        double maxval = Rcpp::max(fv);
        std::vector<double> vc = as<std::vector<double> >(c);
        auto it = std::find_if(vc.begin(), vc.end(), [maxval](double d) { d > maxval; } );
        if (it == vc.begin() || it == vc.end()) Rcpp::stop("This should not happen because candidates should extend beyond patch values"); //TODO
        return Rcpp::wrap(std::vector<double>(it-1,vc.end()));
    }
};

struct AllowDegeneracy {
    Rcpp::NumericVector get_forbidden_values(const DataFrame&) { return Rcpp::NumericVector(); }
    Rcpp::NumericVector filter_borders(const Rcpp::NumericVector& c, const Rcpp::NumericVector&) { return c; }
};



#endif

