#ifndef SCORE_COMPUTER_HPP
#define SCORE_COMPUTER_HPP

#include <Rcpp.h>
#include <tuple> //tie

#include "CalculationTraits.hpp"
#include "util.hpp" //soft_threshold


//Calculation is either Signal or Difference (points to both likelihoods and data structures)
//Score knows how to assemble it into the BIC/CV
template<typename Calculation, typename Score> class ScoreComputer : private Calculation::likelihood_t, private Score {
public:
    ScoreComputer(double tol_val, const typename Calculation::data_t& data, const typename Score::var_t& score_specific) :
       Calculation::likelihood_t(data), Score(score_specific), tol_val_(tol_val), value_(data.get_value()), patchno_(data.get_patchno()) {}
    
    NumericVector evaluate(double LB, double UB) const {
        //compute dof and chi square
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const NumericVector chisq = Calculation::likelihood_t::get_chi_square(LB, UB);
        const int dof = get_dof(LB, UB);
        
        //assemble it into a score
        double score, score_sd;
        std::tie(score, score_sd) = Score::assemble(chisq, dof);
        
        /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
         << " CV= " << score  << " dof= " << dof
         << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
        return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                     _["BIC"]=score, _["BIC.sd"]=score_sd, _["dof"]=dof,
                                     _["UB"]=UB, _["LB"]=LB); // do not report score name for now
    }

private:
    int get_dof(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> value_r = as<std::vector<double> >(value_);
        NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
        IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
        const int dof = unique(selected).size();
        return dof;
    }
    
    double tol_val_;
    NumericVector value_;
    IntegerVector patchno_;

};

#endif
