#ifndef SCORE_COMPUTER_HPP
#define SCORE_COMPUTER_HPP

#include <Rcpp.h>
#include <utility> //pair
#include <tuple> //tie

#include "Traits.hpp"
#include "util.hpp" //soft_threshold
#include "ScoreAssembler.hpp"

//Calculation is either Signal or Difference (points to both likelihoods and data structures)
//Score knows how to assemble it into the BIC/CV
template<typename Calculation, typename Score> class ScoreComputer : private Calculation::likelihood_t, private ScoreAssembler<Score> {
public:
    typedef typename ScoreAssembler<Score>::var_t var_t;
    typedef Rcpp::NumericVector value_t;
    ScoreComputer(double tol_val, const typename Calculation::data_t& data, const var_t& score_specific) :
       Calculation::likelihood_t(data), ScoreAssembler<Score>(score_specific), tol_val_(tol_val), value_(data.get_value()), patchno_(data.get_patchno()) {}
    
    //compute score obtained with these bounds
    Rcpp::NumericVector evaluate(double LB, double UB) const {
        //compute dof and chi square
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const Rcpp::NumericVector chisq = Calculation::likelihood_t::get_chi_square(LB, UB);
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
        const NumericVector dof = NumericVector::get_na();
        /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
         << " CV= " << score  << " dof= " << dof
         << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
        return Rcpp::NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                           _["BIC"]=score, _["BIC.sd"]=score_sd, _["dof"]=dof,
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
