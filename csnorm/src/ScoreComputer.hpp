#ifndef SCORE_COMPUTER_HPP
#define SCORE_COMPUTER_HPP

#include <Rcpp.h>
#include <utility> //pair
#include <tuple> //tie

#include "Traits.hpp"
#include "util.hpp" //soft_threshold
#include "ScoreAssembler.hpp"
#include "Likelihoods.hpp"

//Calculation is either Signal or Difference (points to both likelihoods and data structures)
//Score knows how to assemble it into the BIC/CV
template<typename Calculation, typename Score> class ScoreComputer : private Likelihood<Calculation>, private ScoreAssembler<Score> {
public:
    typedef typename Calculation::binned_t binned_t;
    typedef typename Likelihood<Calculation>::var_t likelihood_var_t;
    typedef typename ScoreAssembler<Score>::var_t assembler_var_t;
    typedef Rcpp::NumericVector value_t;
    ScoreComputer(double tol_val, const binned_t& data, const likelihood_var_t& likelihood_var, const assembler_var_t& assembler_var) :
       Likelihood<Calculation>(data, likelihood_var), ScoreAssembler<Score>(assembler_var), tol_val_(tol_val), value_(data.get_beta()),
       patchno_(data.get_patchno()) {}
    
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
    
private:
    int get_dof(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> value_r = as<std::vector<double> >(value_);
        Rcpp::NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
        Rcpp::IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
        const int dof = unique(selected).size();
        return dof;
    }
    
    double tol_val_;
    Rcpp::NumericVector value_;
    Rcpp::IntegerVector patchno_;

};

#endif
