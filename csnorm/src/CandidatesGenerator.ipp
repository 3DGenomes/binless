#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "graph_helpers.hpp"
#include "BoundsComputer.hpp"
#include "BoundsChecker.hpp"
#include "ScoreComputer.hpp"
#include "ScoreOptimizer.hpp"


template<typename Degeneracy, typename Calculation>
Rcpp::NumericVector
CandidatesGenerator<Degeneracy,Calculation>::expand_values(const NumericVector& values) const {
    std::vector<double> ret;
    ret.reserve(values.size());
    ret.push_back(values(0)-1);
    ret.insert(ret.end(), values.begin(), values.end());
    ret.push_back(values(values.size()-1)+1);
    return Rcpp::wrap(ret);
}

template<typename Degeneracy, typename Calculation>
Rcpp::NumericVector
CandidatesGenerator<Degeneracy,Calculation>::candidates_from_borders(const Rcpp::NumericVector& borders) const {
    std::vector<double> ret;
    if (borders.size()>1) {
        ret.reserve(borders.size()-1);
        for (unsigned i=0; i<borders.size()-1; ++i) {
            ret.push_back((borders[i]+borders[i+1])/2);
        }
    }
    return Rcpp::wrap(ret);
}

template<typename Degeneracy, typename Calculation>
Rcpp::NumericVector
CandidatesGenerator<Degeneracy,Calculation>::get_UB_candidates(int nbins, double tol_val, const binned_t& binned) const {
    //get patch values at which the degrees of freedom change
    Rcpp::NumericVector beta = binned.get_beta(); //Take all patches: even if weight is zero, the dof changes
    Rcpp::IntegerVector patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(), binned.get_bin2(), beta);
    Rcpp::NumericVector borders = get_patch_values(beta, patchno);//TODO: should we operate on beta_cv or beta?
    borders = expand_values(borders);
    //filter out patch values we already know cannot be used to generate UB candidates
    borders = CandidatesFilter<Degeneracy>::filter(borders);
    //compute a list of upper bound candidates
    return candidates_from_borders(borders);
}


