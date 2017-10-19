#ifndef DOF_COMPUTER_HPP
#define DOF_COMPUTER_HPP

#include <Rcpp.h>

#include "Traits.hpp"
#include "util.hpp" //soft_threshold
#include "BinnedData.hpp"

//DOFComputer calculates the degrees of freedom given a set of upper and lower bounds
class DOFComputer {
public:
    DOFComputer(double tol_val, const BinnedDataCore& data) : tol_val_(tol_val), beta_(data.get_beta()),
       patchno_(data.get_patchno()) {}
    
    int get_dof(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> beta_r = as<std::vector<double> >(beta_);
        Rcpp::NumericVector soft = wrap(soft_threshold(beta_r, eCprime, lambda1));
        Rcpp::IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
        const int dof = unique(selected).size();
        return dof;
    }
    
private:
    double tol_val_;
    Rcpp::NumericVector beta_;
    Rcpp::IntegerVector patchno_;

};

#endif
