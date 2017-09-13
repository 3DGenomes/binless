#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "graph_helpers.hpp"
#include "BoundsOptimizer.hpp"
#include "ScoreComputer.hpp"
#include "ScoreOptimizer.hpp"

template<typename Score, typename Offset, typename Sign, typename Degeneracy, typename Calculation, typename GaussianEstimator>
typename SparsityEstimator<Score,Offset,Sign,Degeneracy,Calculation,GaussianEstimator>::score_t
SparsityEstimator<Score,Offset,Sign,Degeneracy,Calculation,GaussianEstimator>::optimize() const {
    //iterate over all candidates and evaluate the score
    std::vector<score_t> values;
    for (double c : CandidatesGenerator<Offset,Sign,Degeneracy>::get_UB_candidates()) {
        //optimize bounds
        bounds_t bounds = BoundsOptimizer<Offset,Sign>::optimize_bounds(c);
        score_t val = ScoreComputer<Calculation,Score,GaussianEstimator>::evaluate(bounds);
        Rcpp::Rcout << "Candidate " << c << " optimized at LB= " << bounds.first << " UB= " << bounds.second
                    << " has score= " << as<double>(val["BIC"]) << " and sd= " << as<double>(val["BIC.sd"]) << "\n";
        values.push_back(val);
    }
    //pick the best one based on the scoring criterion
    score_t best = ScoreOptimizer<Score>::optimize(values); //ScoreOptimizer
    Rcpp::Rcout << "Best is LB= " << as<double>(best["LB"]) << " UB= " << as<double>(best["UB"])
                << " has score= " << as<double>(best["BIC"]) << " and sd= " << as<double>(best["BIC.sd"]) << "\n";
    return best;
}
