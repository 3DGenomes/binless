#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "compute_CV.hpp"
#include "util.hpp"


Rcpp::NumericVector compute_CV_diff::evaluate(double LB, double UB) const {
    //compute dof and CV
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    std::vector<double> value_r = as<std::vector<double> >(value_);
    std::vector<double> soft_r = soft_threshold(value_r, eCprime, lambda1);
    NumericVector soft = wrap(soft_r);
    std::vector<double> valuehat_r = as<std::vector<double> >(valuehat_);
    std::vector<double> valuehat_ref_r = as<std::vector<double> >(valuehat_ref_);
    std::vector<double> var, var_ref;
    var.reserve(valuehat_r.size());
    var_ref.reserve(valuehat_r.size());
    for (int i=0; i<valuehat_r.size(); ++i) {
        var.push_back(1/weight_(i));
        var_ref.push_back(1/weight_ref_(i));
    }
    NumericVector phi_ref = wrap(compute_phi_ref(soft_r, valuehat_r, var, valuehat_ref_r, var_ref));
    
    IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
    const int dof = unique(selected).size();
    
    const NumericVector indiv_CV = weight_ * SQUARE(valuehat_ - (soft + eCprime + phi_ref))
    + weight_ref_ * SQUARE(valuehat_ref_ - phi_ref);
    NumericVector groupwise_CV, groupwise_weights;
    const int ngroups=2;
    for (int i=0; i<ngroups; ++i) {
        groupwise_CV.push_back(sum(as<NumericVector>(indiv_CV[cv_grp_==i])));
        groupwise_weights.push_back(sum(cv_grp_==i));
    }
    const double CV = sum(groupwise_weights*groupwise_CV)/sum(groupwise_weights);
    const double CV_sd = std::sqrt(sum(groupwise_weights*SQUARE(groupwise_CV))/sum(groupwise_weights) - SQUARE(CV));
    
    /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
     << " CV= " << CV  << " dof= " << dof
     << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
    return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                 _["BIC"]=CV, _["BIC.sd"]=CV_sd, _["dof"]=dof,
                                 _["UB"]=UB, _["LB"]=-UB);
}




Rcpp::NumericVector compute_CV_signal::evaluate(double LB, double UB) const {
    //compute dof and CV
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    std::vector<double> value_r = as<std::vector<double> >(value_);
    NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
    IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
    const int dof = unique(selected).size();
    
    const NumericVector indiv_CV = weight_ * SQUARE(valuehat_ - (soft + eCprime));
    NumericVector groupwise_CV, groupwise_weights;
    const int ngroups=2;
    for (int i=0; i<ngroups; ++i) {
        groupwise_CV.push_back(sum(as<NumericVector>(indiv_CV[cv_grp_==i])));
        groupwise_weights.push_back(sum(cv_grp_==i));
    }
    const double CV = sum(groupwise_weights*groupwise_CV)/sum(groupwise_weights);
    const double CV_sd = std::sqrt(sum(groupwise_weights*SQUARE(groupwise_CV))/sum(groupwise_weights) - SQUARE(CV));
    
    /*Rcout << " OBJ " << msg << " ok lambda2= " << lambda2_ << " lambda1= " << lambda1
     << " eCprime= " << eCprime << " CV= " << CV  << " dof= " << dof
     << " UB= " << UB  << " LB= " << LB << std::endl;*/
    return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                 _["BIC"]=CV, _["BIC.sd"]=CV_sd, _["dof"]=dof, _["UB"]=UB, _["LB"]=LB);
}





