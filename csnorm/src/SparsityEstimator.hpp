#ifndef SPARSITY_ESTIMATOR_HPP
#define SPARSITY_ESTIMATOR_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "graph_helpers.hpp"
#include "Traits.hpp"
#include "BoundsComputer.hpp"
#include "BoundsChecker.hpp"
#include "ScoreComputer.hpp"
#include "ScoreOptimizer.hpp"
#include "Degeneracy.hpp"

// A class that estimates lambda1 (and an offset if needed) based on various criteria
/* We define and use upper bound UB = offset + lambda1 and lower bound LB = offset - lambda1
 * BIC and CV have the same quadratic functional form when no patch is added/removed.
 * When UB and LB vary between two patches (or outside any patch), the quadratic has only one minimum
 * In this setting, when UB is chosen, positivity/offset constraints imply the correct LB
 * We thus generate a list of upper bound candidates that get optimized between two patches
 * and evaluate the BIC or the CV at these optimized values
 */
template<typename Calculation, //whether it's a Signal or a Difference calculation
         typename Score, //to choose between BIC, CV or CVkSD
         typename Offset, //whether offset should be held at zero (ZeroOffset) or estimated (EstimatedOffset)
         typename Sign, //whether there is no constraint on the sign of the estimate (AnySign) or whether it must be positive (PositiveSign)
         typename Degeneracy = ForbidDegeneracy> //whether to forbid certain values in order to avoid degeneracies in the model
class SparsityEstimator : private Sign,
                          private Degeneracy,
                          private BoundsComputer<Offset,Sign>,
                          private BoundsChecker<Sign>,
                          private BoundsChecker<Degeneracy>,
                          private ScoreComputer<Calculation,Score>,
                          private ScoreOptimizer<Score> {
    
    typedef typename ScoreComputer<Calculation,Score>::value_t score_t;
    
public:
    
    SparsityEstimator(int nbins, double tol_val, const typename Calculation::data_t& data, double lambda2,
                      const DataFrame& mat, const typename ScoreComputer<Calculation,Score>::var_t& score_specific) :
     Degeneracy(mat),
     BoundsComputer<Offset,Sign>(data, min(data.get_value())),
     BoundsChecker<Sign>(Rcpp::as<Rcpp::NumericVector>(mat["beta"])),
     BoundsChecker<Degeneracy>(Degeneracy::get_forbidden_values()),
     ScoreComputer<Calculation,Score>(tol_val, data, score_specific),
     lambda2_(lambda2), UBcandidates_(get_UB_candidates(nbins, tol_val, mat)) {}
    
    score_t optimize() const {
        //iterate over all candidates and evaluate the score
        std::vector<score_t> values;
        for (double c : UBcandidates_) {
            //optimize bounds
            bounds_t bounds = BoundsComputer<Offset,Sign>::optimize_bounds(c);
            //check if they pass the constraints
            if (BoundsChecker<Sign>::is_valid(bounds) && BoundsChecker<Degeneracy>::is_valid(bounds)) {
                //evaluate the score and store it
                score_t val = ScoreComputer<Calculation,Score>::evaluate(bounds);
                values.push_back(val);
            } else {
                score_t val = ScoreComputer<Calculation,Score>::invalidate(bounds);
                values.push_back(val);
            }
        }
        //pick the best one based on the scoring criterion
        score_t best = ScoreOptimizer<Score>::optimize(values); //ScoreOptimizer
        return best;
    }
    
    
private:
    
    Rcpp::NumericVector expand_values(const NumericVector& values) const {
        std::vector<double> ret;
        ret.reserve(values.size());
        ret.push_back(values(0)-1);
        ret.insert(ret.end(), values.begin(), values.end());
        ret.push_back(values(values.size()-1)+1);
        return Rcpp::wrap(ret);
    }
    
    //propose a set of candidates for an upper bound based on a list of patch borders
    Rcpp::NumericVector candidates_from_borders(const Rcpp::NumericVector& borders) const {
        std::vector<double> ret;
        if (borders.size()>1) {
            ret.reserve(borders.size()-1);
            for (unsigned i=0; i<borders.size()-1; ++i) {
                ret.push_back((borders[i]+borders[i+1])/2);
            }
        }
        return Rcpp::wrap(ret);
    }
    
    Rcpp::NumericVector get_UB_candidates(int nbins, double tol_val, const DataFrame& mat) {
        //get patch values at which the degrees of freedom change
        Rcpp::NumericVector beta = mat["beta"]; //Take all patches: even if weight is zero, the dof changes
        Rcpp::IntegerVector patchno = get_patch_numbers(nbins, mat, tol_val);
        Rcpp::NumericVector borders = get_patch_values(beta, patchno);//TODO: should we operate on beta_cv or beta?
        borders = expand_values(borders);
        //filter out patch values we already know cannot be used to generate UB candidates
        borders = Degeneracy::filter_borders(borders);
        //compute a list of upper bound candidates
        return candidates_from_borders(borders);
    }
    
    const double lambda2_;
    const Rcpp::NumericVector forbidden_values_, UBcandidates_;

};


#endif

