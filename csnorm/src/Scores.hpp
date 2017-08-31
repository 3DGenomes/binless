#ifndef SCORES_HPP
#define SCORES_HPP

#include <Rcpp.h>
#include <utility> //pair

#include "util.hpp" //SQUARE

class CVScore {
public:
    typedef Rcpp::IntegerVector var_t; //type of cv_grp
    typedef std::pair<double,double> value_t; //score return type
    const std::string score_name = "CV";
    CVScore(const var_t& cv_grp) : cv_grp_(cv_grp) {}
    
    value_t assemble(const Rcpp::NumericVector& chisq, double) const {
        Rcpp::NumericVector groupwise_CV, groupwise_weights;
        const int ngroups=2;
        for (int i=0; i<ngroups; ++i) {
            groupwise_CV.push_back(sum(as<NumericVector>(chisq[cv_grp_==i])));
            groupwise_weights.push_back(sum(cv_grp_==i));
        }
        const double CV = sum(groupwise_weights*groupwise_CV)/sum(groupwise_weights);
        const double CV_sd = std::sqrt(sum(groupwise_weights*SQUARE(groupwise_CV))/sum(groupwise_weights) - SQUARE(CV));
        return value_t(CV,CV_sd);
    }
    
private:
    var_t cv_grp_;
};


class BICScore {
public:
    typedef Rcpp::NumericVector var_t; //type of ncounts
    typedef std::pair<double,double> value_t; //score return type
    const std::string score_name = "BIC";
    BICScore(const var_t& ncounts) : lsnc_(log(Rcpp::sum(ncounts))) {}
    
    value_t assemble(const Rcpp::NumericVector& chisq, double dof) const {
        const double BIC = sum(chisq)+ lsnc_*dof;
        const double BIC_sd = -1;
        return value_t(BIC,BIC_sd);
    }
    
private:
    double lsnc_;
};


#endif

