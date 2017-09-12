#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "graph_helpers.hpp"
#include "BoundsComputer.hpp"
#include "BoundsChecker.hpp"
#include "ScoreComputer.hpp"
#include "ScoreOptimizer.hpp"

template<typename Score, typename Offset, typename Sign, typename Degeneracy, typename Calculation, typename GaussianEstimator>
typename SparsityEstimator<Score,Offset,Sign,Degeneracy,Calculation,GaussianEstimator>::score_t
SparsityEstimator<Score,Offset,Sign,Degeneracy,Calculation,GaussianEstimator>::optimize() const {
    //iterate over all candidates and evaluate the score
    std::vector<score_t> values;
    for (double c : UBcandidates_) {
        //optimize bounds
        bounds_t bounds = BoundsComputer<Offset,Sign>::optimize_bounds(c);
        Rcpp::Rcout << "candidate= " << c << " LB= " <<bounds.first << " UB= " << bounds.second;
        //check if they pass the constraints
        if (BoundsChecker<Sign>::is_valid(bounds) && BoundsChecker<Degeneracy>::is_valid(bounds)) {
            //evaluate the score and store it
            score_t val = ScoreComputer<Calculation,Score,GaussianEstimator>::evaluate(bounds);
            Rcpp::Rcout << " score= " << as<double>(val["BIC"]) << " valid\n";
            values.push_back(val);
        } else {
            score_t val = ScoreComputer<Calculation,Score,GaussianEstimator>::invalidate(bounds);
            Rcpp::Rcout << " score= " << as<double>(val["BIC"]) << " invalid sign= " << BoundsChecker<Sign>::is_valid(bounds)
                                      << " degeneracy= " << BoundsChecker<Degeneracy>::is_valid(bounds) << "\n";
            values.push_back(val);
        }
    }
    //pick the best one based on the scoring criterion
    score_t best = ScoreOptimizer<Score>::optimize(values); //ScoreOptimizer
    Rcpp::Rcout << "best score UB= " << as<double>(best["UB"]) << " score= " << as<double>(best["BIC"]) << "\n";
    return best;
}
