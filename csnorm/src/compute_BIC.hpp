#ifndef COMPUTE_BIC_HPP
#define COMPUTE_BIC_HPP

#include "util.hpp"


struct compute_BIC_signal {
    compute_BIC_signal(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
                       const IntegerVector& patchno, const NumericVector& ncounts) :
    tol_val_(tol_val), lsnc_(log(sum(ncounts))), value_(value), weight_(weight), valuehat_(valuehat), patchno_(patchno) {}
    
    NumericVector evaluate(double LB, double UB) const {
        //compute dof and BIC
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> value_r = as<std::vector<double> >(value_);
        NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
        IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
        const int dof = unique(selected).size();
        const double BIC = sum(weight_ * SQUARE(valuehat_ - (soft + eCprime))) + lsnc_*dof;
        /*Rcout << " OBJ " << msg << " ok lambda2= " << lambda2_ << " lambda1= " << lambda1
         << " eCprime= " << eCprime << " BIC= " << BIC  << " dof= " << dof
         << " UB= " << UB  << " LB= " << LB << std::endl;*/
        return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1, _["BIC"]=BIC, _["BIC.sd"]=-1,
                                     _["dof"]=dof, _["UB"]=UB, _["LB"]=LB);
    }

private:
    double tol_val_, lsnc_;
    NumericVector value_, weight_, valuehat_;
    IntegerVector patchno_;
};


#endif