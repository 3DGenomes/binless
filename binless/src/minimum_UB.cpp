#include <Rcpp.h>
#include <algorithm>

#include "minimum_UB.hpp"

template<> double get_minimum_UB_offset_sign<ZeroOffset,PositiveSign>(const BinnedDataCore& binned) {
    return std::max( - (double)(Rcpp::min(binned.get_beta())), 0.);
}
template<> double get_minimum_UB_offset_sign<ZeroOffset,AnySign>(const BinnedDataCore& binned) {
    return 0;
}
template<> double get_minimum_UB_offset_sign<EstimatedOffset,PositiveSign>(const BinnedDataCore& binned) {
    return Rcpp::min(binned.get_beta());
}


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

template<> double get_minimum_UB_degen<AllowDegeneracy>(const BinnedDataCore&, const Rcpp::NumericVector& borders) { return min(borders)-1; }
template<> double get_minimum_UB_degen<ForbidDegeneracy>(const BinnedDataCore& binned, const Rcpp::NumericVector& borders) {
    return get_max_forbidden_border(binned, borders); //folded or not, depending on borders
}
