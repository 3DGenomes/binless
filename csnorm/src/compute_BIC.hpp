#ifndef COMPUTE_BIC_HPP
#define COMPUTE_BIC_HPP

#include "util.hpp"

//DataModel is a temporary placeholder for the object that holds all the data and whether it is a difference or a signal calculation
template<typename DataModel> class compute_BIC : private DataModel {
public:
    //ugly, but works because of SFINAE
    compute_BIC(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
               const IntegerVector& patchno, const NumericVector& ncounts) :
    DataModel(value, weight, valuehat), tol_val_(tol_val), value_(value), patchno_(patchno),
    lsnc_(log(sum(ncounts))) {}
    compute_BIC(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
               const NumericVector& weight_ref, const NumericVector& valuehat_ref,
               const IntegerVector& patchno, const NumericVector& ncounts) :
    DataModel(value, weight, valuehat, weight_ref, valuehat_ref), tol_val_(tol_val),
    value_(value), patchno_(patchno), lsnc_(log(sum(ncounts))) {}
    
    NumericVector evaluate(double LB, double UB) const {
        //compute dof and CV
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        const NumericVector chisq = DataModel::get_chi_square(LB, UB);
        const int dof = get_dof(LB, UB);

        const double BIC = sum(chisq)+ lsnc_*dof;
        const double BIC_sd = -1;

        /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
         << " CV= " << CV  << " dof= " << dof
         << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
        return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                     _["BIC"]=BIC, _["BIC.sd"]=BIC_sd, _["dof"]=dof,
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
    IntegerVector patchno_;
    double lsnc_;
    
};



#endif
