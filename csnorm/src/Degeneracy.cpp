#include <Rcpp.h>
#include <vector>
#include <algorithm>

#include "Degeneracy.hpp"


Rcpp::NumericVector ForbidDegeneracy::filter_borders(const Rcpp::NumericVector& c) const {
    double maxval = Rcpp::max(forbidden_values_);
    std::vector<double> vc = Rcpp::as<std::vector<double> >(c);
    auto it = std::find_if(vc.begin(), vc.end(), [maxval](double d) { return d > maxval; } );
    if (it == vc.begin() || it == vc.end()) Rcpp::stop("This should not happen because candidates should extend beyond patch values"); //TODO
    return Rcpp::wrap(std::vector<double>(it-1,vc.end()));
}
    
Rcpp::NumericVector ForbidDegeneracy::compute_forbidden_values(const Rcpp::DataFrame& mat) const {
    Rcpp::NumericVector beta = mat["beta"];
    Rcpp::IntegerVector diag_grp = mat["diag.grp"];
    return get_minimum_diagonal_values(beta, diag_grp); //TODO: beta or abs(beta) or a combination?
}

Rcpp::NumericVector ForbidDegeneracy::get_minimum_diagonal_values(const Rcpp::NumericVector& value, const Rcpp::IntegerVector& diag_grp) const {
    int ndiags = max(diag_grp)+1;
    Rcpp::LogicalVector diagtouch(ndiags, false);
    Rcpp::NumericVector diagvals(ndiags, max(value)); //diag_idx starts at 0
    for (int i=0; i<diag_grp.size(); ++i) {
        diagvals(diag_grp(i)) = std::min(diagvals(diag_grp(i)),value(i));
        diagtouch(diag_grp(i)) = true;
    }
    diagvals = diagvals[diagtouch];
    diagvals = unique(diagvals);
    std::sort(diagvals.begin(), diagvals.end());
    return diagvals;
}

