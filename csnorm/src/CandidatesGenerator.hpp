#ifndef CANDIDATES_GENERATOR_HPP
#define CANDIDATES_GENERATOR_HPP

#include <Rcpp.h>
#include <vector>

#include "CandidatesFilter.hpp"
#include "BinnedData.hpp"

template<typename Degeneracy, //whether to forbid certain values in order to avoid degeneracies in the model
         typename Calculation> //whether it's a Signal or a Difference calculation
class CandidatesGenerator : private CandidatesFilter<Degeneracy> {
    
    typedef BinnedData<Calculation> binned_t;
    
public:
    
    CandidatesGenerator(const binned_t& binned) : CandidatesFilter<Degeneracy>(binned) {}
    
    Rcpp::NumericVector get_UB_candidates(int nbins, double tol_val, const binned_t& binned) const;
};

Rcpp::NumericVector expand_values(const Rcpp::NumericVector& values);

//propose a set of candidates for an upper bound based on a list of patch borders
Rcpp::NumericVector candidates_from_borders(const Rcpp::NumericVector& borders);

Rcpp::NumericVector get_patch_values(const Rcpp::NumericVector& value, const Rcpp::IntegerVector& patchno);


#include "CandidatesGenerator.ipp"


#endif

