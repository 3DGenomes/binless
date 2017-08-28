#ifndef COMPUTE_CV_HPP
#define COMPUTE_CV_HPP

#include "util.hpp"

//DataModel is a temporary placeholder for the object that holds all the data and whether it is a difference or a signal calculation
template<typename DataModel> class compute_CV : private DataModel {
public:
    //ugly, but works because of SFINAE
    compute_CV(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
               const IntegerVector& patchno, const IntegerVector& cv_grp) :
       DataModel(value, weight, valuehat), tol_val_(tol_val), value_(value), patchno_(patchno), cv_grp_(cv_grp) {}
    compute_CV(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
               const NumericVector& weight_ref, const NumericVector& valuehat_ref,
               const IntegerVector& patchno, const IntegerVector& cv_grp) :
       DataModel(value, weight, valuehat, weight_ref, valuehat_ref), tol_val_(tol_val), value_(value), patchno_(patchno), cv_grp_(cv_grp) {}
    
    NumericVector evaluate(double LB, double UB) const {
        //compute dof and CV
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const NumericVector chisq = DataModel::get_chi_square(LB, UB);
        
        NumericVector groupwise_CV, groupwise_weights;
        const int ngroups=2;
        for (int i=0; i<ngroups; ++i) {
            groupwise_CV.push_back(sum(as<NumericVector>(chisq[cv_grp_==i])));
            groupwise_weights.push_back(sum(cv_grp_==i));
        }
        const double CV = sum(groupwise_weights*groupwise_CV)/sum(groupwise_weights);
        const double CV_sd = std::sqrt(sum(groupwise_weights*SQUARE(groupwise_CV))/sum(groupwise_weights) - SQUARE(CV));
        
        const int dof = get_dof(LB, UB);
        /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
         << " CV= " << CV  << " dof= " << dof
         << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
        return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                     _["BIC"]=CV, _["BIC.sd"]=CV_sd, _["dof"]=dof,
                                     _["UB"]=UB, _["LB"]=-UB);
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
    IntegerVector patchno_, cv_grp_;

};

#endif