#include <Rcpp.h>
#include <vector>
#include <algorithm>

#include "CandidatesGenerator.hpp"

Rcpp::NumericVector expand_values(const NumericVector& values) {
    std::vector<double> ret;
    ret.reserve(values.size());
    ret.push_back(values(0)-1);
    ret.insert(ret.end(), values.begin(), values.end());
    ret.push_back(values(values.size()-1)+1);
    return Rcpp::wrap(ret);
}

Rcpp::NumericVector candidates_from_borders(const Rcpp::NumericVector& borders) {
    std::vector<double> ret;
    if (borders.size()>1) {
        ret.reserve(borders.size()-1);
        for (unsigned i=0; i<borders.size()-1; ++i) {
            ret.push_back((borders[i]+borders[i+1])/2);
        }
    }
    return Rcpp::wrap(ret);
}

//return ascending list of the values of all detected patches in the matrix
Rcpp::NumericVector get_patch_values(const NumericVector& value, const IntegerVector& patchno) {
    int npatches = max(patchno)+1;
    NumericVector unique_values(npatches); //patchno starts at 0
    for (int i=0; i<patchno.size(); ++i) unique_values(patchno(i)) = value(i);
    std::sort(unique_values.begin(), unique_values.end());
    return unique_values;
}

