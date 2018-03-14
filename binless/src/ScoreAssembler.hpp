#ifndef SCORE_ASSEMBLER_HPP
#define SCORE_ASSEMBLER_HPP

#include <Rcpp.h>
#include <utility> //pair

#include "util.hpp" //SQUARE
#include "Traits.hpp"

//ScoreAssembler receives a chi square and the degrees of freedom
//and assembles them into a score (templated, BIC or CVkSD)
template<typename Score> class ScoreAssembler {};

//CV with or without kSD
template<unsigned kSD> class ScoreAssembler<CVkSD<kSD> > {
public:
    typedef Rcpp::IntegerVector var_t; //type of cv_grp
    typedef std::pair<double,double> value_t; //score return type
    const std::string score_name = "CV";
    
    ScoreAssembler(const var_t& cv_grp) : cv_grp_(cv_grp) {}
    
    value_t assemble(const Rcpp::NumericVector& chisq, double) const {
        Rcpp::NumericVector groupwise_CV, groupwise_weights;
        const int ngroups=2;
        for (int i=0; i<ngroups; ++i) {
            double Ni = sum(cv_grp_==i);
            groupwise_weights.push_back(Ni);
            groupwise_CV.push_back(sum(as<NumericVector>(chisq[cv_grp_==i]))/Ni);
        }
        const double N = sum(groupwise_weights);
        const double CV = sum(groupwise_weights*groupwise_CV)/N;
        const double CV_sd = std::sqrt(sum(groupwise_weights*SQUARE(groupwise_CV-CV))/(N-1));
        return value_t(CV,CV_sd);
    }
    
private:
    var_t cv_grp_;
};

//BIC
template<> class ScoreAssembler<BIC> {
public:
    typedef Rcpp::NumericVector var_t; //type of nobs
    typedef std::pair<double,double> value_t; //score return type
    const std::string score_name = "BIC";
    
    ScoreAssembler(const var_t& nobs) : lsnc_(log(Rcpp::sum(nobs))) {}
    
    value_t assemble(const Rcpp::NumericVector& chisq, double dof) const {
        const double BIC = sum(chisq)+ lsnc_*dof;
        const double BIC_sd = -1;
        return value_t(BIC,BIC_sd);
    }
    
private:
    double lsnc_;
};


#endif

