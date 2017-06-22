#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#define SQUARE(x) ((x)*(x))

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1);

NumericVector get_patch_values(NumericVector value, IntegerVector patchno);

NumericVector get_minimum_diagonal_values(NumericVector value,
        IntegerVector diag_idx);

NumericVector get_constant_diagonal_values(NumericVector value,
        IntegerVector diag_idx, double tol_val);

#endif

