#include <Rcpp.h>
#include <vector>
#include <algorithm>

#include "graph_helpers.hpp" //get_patch_numbers
#include "candidate_generation.hpp"

//return ascending list of the values of all detected patches in the matrix
Rcpp::NumericVector get_patch_values(const Rcpp::NumericVector& value, const Rcpp::IntegerVector& patchno) {
    int npatches = max(patchno)+1;
    Rcpp::NumericVector unique_values(npatches); //patchno starts at 0
    for (int i=0; i<patchno.size(); ++i) unique_values(patchno(i)) = value(i);
    std::sort(unique_values.begin(), unique_values.end());
    return unique_values;
}

Rcpp::NumericVector borders_from_values(int nbins, double tol_val, const BinnedDataCore& binned, const Rcpp::NumericVector& values) {
    Rcpp::IntegerVector patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(), binned.get_bin2(), values);
    Rcpp::NumericVector patchvals = get_patch_values(values, patchno);
    return patchvals;
}

Rcpp::NumericVector cleanup_borders(const Rcpp::NumericVector& borders, double minUB) {
    Rcpp::NumericVector newborders = borders[borders>=minUB];
    std::vector<double> ret;
    ret.reserve(newborders.size());
    ret.push_back(minUB);
    ret.insert(ret.end(), newborders.begin(), newborders.end());
    if (newborders.size()==0) {
        ret.push_back(minUB+1);
    } else {
        ret.push_back(newborders(newborders.size()-1)+1);
    }
    return Rcpp::wrap(ret);
}

Rcpp::NumericVector candidates_from_borders(const Rcpp::NumericVector& borders) {
    std::vector<double> ret;
    ret.reserve(borders.size()-1); //borders is at least size 2, see cleanup_borders
    for (unsigned i=0; i<borders.size()-1; ++i) {
        ret.push_back((borders[i]+borders[i+1])/2);
    }
    return Rcpp::wrap(ret);
}

