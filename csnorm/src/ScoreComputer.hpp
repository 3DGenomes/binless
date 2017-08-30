#ifndef SCORE_COMPUTER_HPP
#define SCORE_COMPUTER_HPP

#include <Rcpp.h>
#include <tuple> //tie

#include "util.hpp" //soft_threshold


//DataLikelihood knows how to compute the chi square, Score knows how to assemble it into the BIC/CV
template<typename DataLikelihood, typename Score> class ScoreComputer : private DataLikelihood, private Score {
public:
    //ugly, but works because of SFINAE
    ScoreComputer(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
                  const IntegerVector& patchno, const typename Score::var_t& score_specific) :
       DataLikelihood(value, weight, valuehat), Score(score_specific), tol_val_(tol_val), value_(value), patchno_(patchno) {}
    
    template<typename Data>
    ScoreComputer(double tol_val, const Data& data, const typename Score::var_t& score_specific) :
       DataLikelihood(data), Score(score_specific), tol_val_(tol_val), value_(data.get_value()), patchno_(data.get_patchno()) {}
    
    ScoreComputer(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
               const NumericVector& weight_ref, const NumericVector& valuehat_ref,
                  const IntegerVector& patchno, const typename Score::var_t& score_specific) :
       DataLikelihood(value, weight, valuehat, weight_ref, valuehat_ref), Score(score_specific), tol_val_(tol_val), value_(value), patchno_(patchno) {}
    
    NumericVector evaluate(double LB, double UB) const {
        //compute dof and chi square
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const NumericVector chisq = DataLikelihood::get_chi_square(LB, UB);
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
