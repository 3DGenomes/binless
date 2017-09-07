#ifndef SCORE_COMPUTER_HPP
#define SCORE_COMPUTER_HPP

#include <Rcpp.h>
#include <utility> //pair
#include <tuple> //tie

#include "Traits.hpp"
#include "ScoreAssembler.hpp"
#include "Likelihoods.hpp"
#include "DOFComputer.hpp"
#include "BinnedData.hpp"

//ScoreComputer takes data and evaluates the score at a given set of upper and lower bounds
//Calculation is either Signal or Difference (points to both likelihoods and data structures)
//Score knows how to assemble it into the BIC/CV
template<typename Calculation, typename Score>
class ScoreComputer : private Likelihood<Calculation>,
                      private ScoreAssembler<Score>,
                      private DOFComputer {
public:
    typedef BinnedData<Calculation> binned_t;
    typedef typename Likelihood<Calculation>::var_t likelihood_var_t;
    typedef typename ScoreAssembler<Score>::var_t assembler_var_t;
    typedef Rcpp::NumericVector value_t;
    ScoreComputer(double tol_val, const binned_t& data, const likelihood_var_t& likelihood_var, const assembler_var_t& assembler_var);
    
    //compute score obtained with these bounds
    Rcpp::NumericVector evaluate(double LB, double UB) const;
    Rcpp::NumericVector evaluate(std::pair<double,double> bounds) const { return evaluate(bounds.first, bounds.second); }

                          
    //return highest possible score, invalidating these bounds
    Rcpp::NumericVector invalidate(double LB, double UB) const;
    Rcpp::NumericVector invalidate(std::pair<double,double> bounds) const { return invalidate(bounds.first, bounds.second); }
};

#include "ScoreComputer.ipp"

#endif
