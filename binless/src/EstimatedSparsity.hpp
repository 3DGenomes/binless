#ifndef ESTIMATED_SPARSITY_HPP
#define ESTIMATED_SPARSITY_HPP

#include <Rcpp.h>
#include <vector>

#include "BoundsOptimizer.hpp"
#include "ScoreComputer.hpp"
#include "ScoreOptimizer.hpp"
#include "util.hpp"
#include "BinnedData.hpp"
#include "candidate_generation.hpp"

// A class that estimates lambda1 (and an offset if needed) based on various criteria
/* We define and use upper bound UB = offset + lambda1 and lower bound LB = offset - lambda1
 * BIC and CV have the same quadratic functional form when no patch is added/removed.
 * When UB and LB vary between two patches (or outside any patch), the quadratic has only one minimum
 * In this setting, when UB is chosen, positivity/offset constraints imply the correct LB
 * We thus generate a list of upper bound candidates that get optimized between two patches
 * and evaluate the BIC or the CV at these optimized values
 */
template<typename Score, //to choose between BIC, CV or CVkSD
         typename Offset, //whether offset should be held at zero (ZeroOffset) or estimated (EstimatedOffset)
         typename Sign, //whether there is no constraint on the sign of the estimate (AnySign) or whether it must be positive (PositiveSign)
         typename Degeneracy, //whether to forbid certain values in order to avoid degeneracies in the model
         typename Calculation, //whether it's a Signal or a Difference calculation
         typename GaussianEstimator> //the type of the fused lasso object used for initialization
class EstimatedSparsity : private BoundsOptimizer<Offset,Sign>,
                          private ScoreComputer<Calculation,Score,GaussianEstimator>,
                          private ScoreOptimizer<Score> {
    
    typedef typename ScoreComputer<Calculation,Score,GaussianEstimator>::value_t score_t;
    typedef BinnedData<Calculation> binned_t;
    
public:
    
    EstimatedSparsity(int nbins, double tol_val, const binned_t& binned, double lambda2,
                      GaussianEstimator& gauss) :
     BoundsOptimizer<Offset,Sign>(binned),
     ScoreComputer<Calculation,Score,GaussianEstimator>(tol_val, binned, gauss, lambda2),
     ScoreOptimizer<Score>(),
     UBcandidates_(get_UB_candidates<Offset,Sign,Degeneracy>(nbins, tol_val, binned)),
     lambda2_(lambda2), lambda1_min_(tol_val) {}
    
     score_t optimize() const;
    
private:
    const Rcpp::NumericVector UBcandidates_;
    const double lambda2_, lambda1_min_;
};

//named constructor
template<typename Score, typename Offset, typename Sign, typename Degeneracy, typename Calculation, typename GaussianEstimator>
EstimatedSparsity<Score, Offset, Sign, Degeneracy, Calculation, GaussianEstimator>
make_EstimatedSparsity(int nbins, double tol_val, const BinnedData<Calculation>& binned, double lambda2,
                       GaussianEstimator& gauss) {
    return EstimatedSparsity<Score, Offset, Sign, Degeneracy, Calculation, GaussianEstimator>(nbins, tol_val, binned, lambda2, gauss);
}


#include "EstimatedSparsity.ipp"


#endif

