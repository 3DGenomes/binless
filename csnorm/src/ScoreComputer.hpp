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
    ScoreComputer(double tol_val, const binned_t& data, const likelihood_var_t& likelihood_var, const assembler_var_t& assembler_var) :
       Likelihood<Calculation>(data, likelihood_var), ScoreAssembler<Score>(assembler_var), DOFComputer(tol_val, data) {}
    
    //compute score obtained with these bounds
    Rcpp::NumericVector evaluate(double LB, double UB) const {
        //compute dof and chi square
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const Rcpp::NumericVector chisq = Likelihood<Calculation>::get_chi_square(LB, UB);
        const int dof = get_dof(LB, UB);
        
        //assemble it into a score
        double score, score_sd;
        std::tie(score, score_sd) = ScoreAssembler<Score>::assemble(chisq, dof);
        
        /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
         << " CV= " << score  << " dof= " << dof
         << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
        return Rcpp::NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                     _["BIC"]=score, _["BIC.sd"]=score_sd, _["dof"]=dof,
                                     _["UB"]=UB, _["LB"]=LB); // do not report score name for now
    }

    Rcpp::NumericVector evaluate(std::pair<double,double> bounds) const { return evaluate(bounds.first, bounds.second); }
   
    //return highest possible score, invalidating these bounds
    Rcpp::NumericVector invalidate(double LB, double UB) const {
        //compute dof and chi square
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const double score = std::numeric_limits<double>::max();
        const double score_sd = -1;
        /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
         << " CV= " << score  << " dof= NA"
         << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
        return Rcpp::NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                           _["BIC"]=score, _["BIC.sd"]=score_sd, _["dof"]=NumericVector::get_na(),
                                           _["UB"]=UB, _["LB"]=LB); // do not report score name for now
    }
    
    Rcpp::NumericVector invalidate(std::pair<double,double> bounds) const { return invalidate(bounds.first, bounds.second); }
};

#endif
