#ifndef CANDIDATES_GENERATOR_HPP
#define CANDIDATES_GENERATOR_HPP

#include <Rcpp.h>
#include <vector>
#include <algorithm>

#include "CandidatesValues.hpp"
#include "CandidatesUBMinDegen.hpp"
#include "CandidatesUBMinOsi.hpp"
#include "BinnedData.hpp"



Rcpp::NumericVector get_patch_values(const Rcpp::NumericVector& value, const Rcpp::IntegerVector& patchno);
Rcpp::NumericVector borders_from_values(int nbins, double tol_val, const BinnedDataCore& binned, const Rcpp::NumericVector& values);

Rcpp::NumericVector cleanup_borders(const Rcpp::NumericVector& borders, double minUB);

//propose a set of candidates for an upper bound based on a list of patch borders
Rcpp::NumericVector candidates_from_borders(const Rcpp::NumericVector& borders);



template<typename Offset, //whether offset should be held at zero (ZeroOffset) or estimated (EstimatedOffset)
         typename Sign, //whether there is no constraint on the sign of the estimate (AnySign) or whether it must be positive (PositiveSign)
         typename Degeneracy> //whether to forbid certain values (ForbidDegeneracy) in order to avoid degeneracies in the model, or not (AllowDegeneracy)
class CandidatesGenerator : private CandidatesValues<Offset,Sign>,
                            private CandidatesUBMinDegen<Degeneracy>,
                            private CandidatesUBMinOsi<Offset,Sign> {
    
public:
    
    CandidatesGenerator(int nbins, double tol_val, const BinnedDataCore& binned)
      : CandidatesValues<Offset,Sign>(binned),
        CandidatesUBMinDegen<Degeneracy>(binned),
        CandidatesUBMinOsi<Offset,Sign>(binned),
        UBcandidates_(generate(nbins, tol_val, binned)) {}

    Rcpp::NumericVector get_UB_candidates() const { return UBcandidates_; }

private:
    Rcpp::NumericVector generate(int nbins, double tol_val, const BinnedDataCore& binned) const {
        //get UB borders: values at which dof changes (including zero weight)
        Rcpp::Rcout << "beta min=" << min(binned.get_beta()) << " max=" << max(binned.get_beta()) << " N=" << binned.get_beta().size() << "\n";
        Rcpp::NumericVector values = CandidatesValues<Offset,Sign>::get();
        Rcpp::Rcout << "values min=" << min(values) << " max=" << max(values) << " N=" << values.size() << "\n";
        Rcpp::NumericVector borders = borders_from_values(nbins,tol_val, binned, values);
        Rcpp::Rcout << "borders min=" << min(borders) << " max=" << max(borders) << " N=" << borders.size() << "\n";
        //get minimum admissible UB
        double minUBd = CandidatesUBMinDegen<Degeneracy>::get(values);
        double minUBo = CandidatesUBMinOsi<Offset,Sign>::get();
        double minUB = std::max(minUBd,minUBo);
        Rcout<< "minUB = max("<<minUBd<<", " << minUBo << ")\n";
        //cleanup the borders list, adding minUB as the lowest possible, expanding at the end
        borders = cleanup_borders(borders, minUB);
        Rcpp::Rcout << "cleanup min=" << min(borders) << " max=" << max(borders) << " N=" << borders.size() << "\n";
        //now generate candidates from these borders
        //pick center of each interval to make it an UB candidate
        Rcpp::NumericVector UBcandidates = candidates_from_borders(borders);
        Rcpp::Rcout << "candidates min=" << min(UBcandidates) << " max=" << max(UBcandidates) << " N=" << UBcandidates.size() << "\n";
        return UBcandidates;
    }
    
    const Rcpp::NumericVector UBcandidates_;
};

#endif

