#ifndef FIXED_SPARSITY_HPP
#define FIXED_SPARSITY_HPP

#include <Rcpp.h>
#include <vector>

#include "BoundsRectifier.hpp"
#include "ScoreComputer.hpp"
#include "util.hpp"
#include "BinnedData.hpp"
#include "candidate_generation.hpp"

// A class that proposes lambda1 (and an offset if needed) as close as possible to the one provided
template<typename Score, //to choose between BIC, CV or CVkSD
         typename Offset, //whether offset should be held at zero (ZeroOffset) or estimated (EstimatedOffset)
         typename Sign, //whether there is no constraint on the sign of the estimate (AnySign) or whether it must be positive (PositiveSign)
         typename Degeneracy, //whether to forbid certain values in order to avoid degeneracies in the model
         typename Calculation, //whether it's a Signal or a Difference calculation
         typename GaussianEstimator> //the type of the fused lasso object used for initialization
class FixedSparsity : private BoundsRectifier<Offset,Sign>,
                      private ScoreComputer<Calculation,Score,GaussianEstimator> {
    
    typedef typename ScoreComputer<Calculation,Score,GaussianEstimator>::value_t score_t;
    typedef BinnedData<Calculation> binned_t;
    
public:
    
    FixedSparsity(int nbins, double tol_val, const binned_t& binned, double lambda2,
                      GaussianEstimator& gauss, double lambda1_candidate) :
     BoundsRectifier<Offset,Sign>(binned),
     ScoreComputer<Calculation,Score,GaussianEstimator>(tol_val, binned, gauss, lambda2),
      lambda1_candidate_(lambda1_candidate),
      lambda2_(lambda2), minUB_(get_minimum_UB<Offset,Sign,Degeneracy>(binned)) {}
    
    score_t optimize() const {
        //convert lambda1 to the closest admissible bounds
        auto bounds = BoundsRectifier<Offset,Sign>::rectify_bounds(lambda1_candidate_, minUB_);
        //score
        score_t val = ScoreComputer<Calculation,Score,GaussianEstimator>::evaluate(bounds);
        /*Rcpp::Rcout << "FixedSparsity: lambda1_candidate= " << lambda1_candidate_ << " minUB_= " << minUB_
            << " lambda1= " << Rcpp::as<double>(val["lambda1"]) <<  " eCprime= " << Rcpp::as<double>(val["eCprime"])
            << " LB= " << bounds.first <<  " UB= " << bounds.second
            << " val= " << Rcpp::as<double>(val["BIC"]) << " dof= " << Rcpp::as<unsigned>(val["dof"]) << "\n";*/
        return val;
    }

    
private:
    const double lambda1_candidate_;
    const double lambda2_;
    const double minUB_;
};

//named constructor
template<typename Score, typename Offset, typename Sign, typename Degeneracy, typename Calculation, typename GaussianEstimator>
FixedSparsity<Score, Offset, Sign, Degeneracy, Calculation, GaussianEstimator>
make_FixedSparsity(int nbins, double tol_val, const BinnedData<Calculation>& binned, double lambda2,
                       GaussianEstimator& gauss, double lambda1_candidate) {
    return FixedSparsity<Score, Offset, Sign, Degeneracy, Calculation, GaussianEstimator>(nbins, tol_val, binned, lambda2, gauss, lambda1_candidate);
}


#endif

