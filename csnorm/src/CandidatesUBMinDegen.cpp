#include <Rcpp.h>
#include <algorithm>

#include "CandidatesUBMinDegen.hpp"


double get_max_forbidden_border(const BinnedDataCore& binned, const Rcpp::NumericVector& beta) {
    Rcpp::IntegerVector diag_grp = binned.get_diag_grp();
    //retrieve minimum values for each counter diagonal group
    int ndiags = max(diag_grp)+1;
    Rcpp::LogicalVector diagtouch(ndiags, false);
    Rcpp::NumericVector diagvals(ndiags, max(beta)); //diag_idx starts at 0
    for (int i=0; i<diag_grp.size(); ++i) {
        diagvals(diag_grp(i)) = std::min(diagvals(diag_grp(i)),beta(i));
        diagtouch(diag_grp(i)) = true;
    }
    diagvals = diagvals[diagtouch];
    diagvals = unique(diagvals);
    return max(diagvals);
}
