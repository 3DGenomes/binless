#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include "util.hpp"

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1) {
    std::vector<double> phi;
    phi.reserve(beta.size());
    for (std::vector<double>::const_iterator it = beta.begin(); it != beta.end();
            ++it) {
        double val = *it - eCprime;
        phi.push_back((val > 0) ? std::max(0., val - lam1) : std::min(0., val + lam1));
    }
    return phi;
}

//return ascending list of the values of all detected patches in the matrix
NumericVector get_patch_values(NumericVector value, IntegerVector patchno) {
    int npatches = max(patchno)+1;
    NumericVector unique_values(npatches); //patchno starts at 0
    for (int i=0; i<patchno.size(); ++i) unique_values(patchno(i)) = value(i);
    std::sort(unique_values.begin(), unique_values.end());
    return unique_values;
}

//return the minimum value encountered along each counter diagonal
NumericVector get_minimum_diagonal_values(NumericVector value,
        IntegerVector diag_idx) {
    int ndiags = max(diag_idx)+1;
    NumericVector diagvals(ndiags, max(value)); //diag_idx starts at 0
    for (int i=0; i<diag_idx.size();
            ++i) diagvals(diag_idx(i)) = std::min(diagvals(diag_idx(i)),value(i));
    return diagvals;
}

//return those values whose patches contain a whole counter diagonal
NumericVector get_constant_diagonal_values(NumericVector value,
        IntegerVector diag_idx, double tol_val) {
    int ndiags = max(diag_idx)+1;
    NumericVector diag_min(ndiags, max(value)), diag_max(ndiags,
            min(value)); //diag_idx starts at 0
    for (int i=0; i<diag_idx.size(); ++i) {
        diag_min(diag_idx(i)) = std::min(diag_min(diag_idx(i)),value(i));
        diag_max(diag_idx(i)) = std::max(diag_max(diag_idx(i)),value(i));
    }
    NumericVector diagvals = (diag_max+diag_min)/2.;
    diagvals = diagvals[(diag_max-diag_min)<=tol_val];
    return diagvals;
}



