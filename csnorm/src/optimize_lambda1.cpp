#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "optimize_lambda1.hpp"
#include "util.hpp"
#include "graph_helpers.hpp" //get_patch_numbers



Rcpp::NumericVector compute_BIC_signal::evaluate(double LB, double UB) const {
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



obj_lambda1_base::obj_lambda1_base(NumericVector value, NumericVector weight, NumericVector valuehat, double minUB) :
    value_(value), absval_(abs(value)), weight_(weight), valuehat_(valuehat), minabsval_(min(absval_)),
    maxabsval_(max(absval_)), minUB_(minUB) {}

obj_lambda1_base::bounds_t obj_lambda1_base::optimize_bounds(double val) const {
    //split data in two groups and determine constraint values
    //Rcout << "  GET at " << val << " maxabsval_= " << maxabsval_ << std::endl;
    LogicalVector grp1 = value_ > val, grp2 = value_ < -val;
    double xk,xkp1;
    if (val>maxabsval_) {
        xk=maxabsval_;
        xkp1=std::numeric_limits<double>::infinity();
    } else if (val<minabsval_) {
        xk=0;
        xkp1=minabsval_;
    } else {
        xk=max(as<NumericVector>(absval_[absval_<=val]));
        xkp1=min(as<NumericVector>(absval_[absval_>val]));
    }
    //Rcout << "  xk= " << xk << " xkp1= " << xkp1 << std::endl;
    //determine unconstrained UB
    double UB;
    bool grp1empty = is_true(all(!grp1)), grp2empty = is_true(all(!grp2));
    //Rcout << "  grp1empty= " << grp1empty << " grp2empty= " << grp2empty << std::endl;
    if (grp1empty && grp2empty) {
        //UB larger than any value
        UB=maxabsval_;
    } else if (grp1empty) {
        //UB larger than any positive value
        NumericVector w2 = weight_[grp2];
        NumericVector betahat2 = valuehat_[grp2];
        NumericVector x2 = value_[grp2];
        UB = -sum(w2*(x2-betahat2))/sum(w2);
    } else if (grp2empty) {
        NumericVector w1 = weight_[grp1];
        NumericVector betahat1 = valuehat_[grp1];
        NumericVector x1 = value_[grp1];
        UB = sum(w1*(x1-betahat1))/sum(w1);
    } else {
        NumericVector w1 = weight_[grp1];
        NumericVector betahat1 = valuehat_[grp1];
        NumericVector x1 = value_[grp1];
        double sw1 = sum(w1);
        double b1 = sum(w1*(x1-betahat1))/sw1;
        NumericVector w2 = weight_[grp2];
        NumericVector betahat2 = valuehat_[grp2];
        NumericVector x2 = value_[grp2];
        double sw2 = sum(w2);
        double b2 = sum(w2*(x2-betahat2))/sw2;
        UB = (sw1*b1-sw2*b2)/(sw1+sw2);
    }
    //apply constraint
    UB = std::min(std::max(UB,xk),xkp1);
    if (minUB_ <= xkp1 && minUB_ > xk) UB=std::max(minUB_,UB);
    //Rcout << "  UB= " << UB << " minUB= " << minUB_ << std::endl;
    return UB;
}

    
    
    

obj_lambda1_BIC::obj_lambda1_BIC(double minUB, double tol_val,
                                 IntegerVector patchno, NumericVector forbidden_vals,
                                 NumericVector value, NumericVector weight, NumericVector valuehat,
                                 NumericVector ncounts) :
  obj_lambda1_base(value, weight, valuehat, minUB),
  compute_BIC_signal(tol_val, value, weight, valuehat, patchno, ncounts),
  minUB_(minUB), tol_val_(tol_val), lsnc_(log(sum(ncounts))), forbidden_vals_(forbidden_vals),
  patchno_(patchno), value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_BIC::operator()(double x) const {
  //return get(std::pow(10,x), "opt")["BIC"];
  return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_BIC::get(double val, std::string msg) const {
  
    double UB = optimize_bounds(val);
    
  //check if forbidden
  double lambda1=UB;
  if ( is_true(any(abs(forbidden_vals_)>lambda1+tol_val_/2)) || (UB < minUB_) ) {
    if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1
                            << " eCprime= 0 BIC= Inf dof= NA"
                            << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
    return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                 _["BIC"]=std::numeric_limits<double>::max(), _["dof"]=NumericVector::get_na(),
                                 _["UB"]=lambda1, _["LB"]=-lambda1);
  }
    return evaluate(-UB, UB);
}




obj_lambda1_CV::obj_lambda1_CV(double minUB, double tol_val,
                         IntegerVector patchno, NumericVector forbidden_vals,
                         NumericVector value, NumericVector weight, NumericVector valuehat,
                         NumericVector ncounts, IntegerVector cv_grp) :
    obj_lambda1_base(value, weight, valuehat, minUB),
    compute_CV_signal(tol_val, value, weight, valuehat, patchno, cv_grp),
    minUB_(minUB), tol_val_(tol_val), lsnc_(log(sum(ncounts))), forbidden_vals_(forbidden_vals),
    patchno_(patchno), value_(value), weight_(weight), valuehat_(valuehat), cv_grp_(cv_grp) {}

double obj_lambda1_CV::operator()(double x) const {
    //return get(std::pow(10,x), "opt")["BIC"];
    return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_CV::get(double val, std::string msg) const {
    
    double UB = optimize_bounds(val);
    
    //check if forbidden
  double lambda1=UB;
  if ( is_true(any(abs(forbidden_vals_)>lambda1+tol_val_/2)) || (UB < minUB_) ) {
      if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1
                              << " eCprime= 0 CV= Inf dof= NA"
                              << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
      return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                   _["BIC"]=std::numeric_limits<double>::max(), _["BIC.sd"]=0,
                                   _["dof"]=NumericVector::get_na(), _["UB"]=lambda1, _["LB"]=-lambda1);
  }
  
  return evaluate(-UB, UB);
}

NumericVector cpp_optimize_lambda1(const DataFrame mat, int nbins,
                                   double tol_val, bool positive,
                                   double lambda1_min, int refine_num) {
    const bool constrained = true;
    //extract vectors
    double lmin = std::max(lambda1_min,tol_val/2);
    NumericVector weight = mat["weight"];
    NumericVector phihat = mat["phihat"];
    NumericVector beta = mat["beta"];
    NumericVector ncounts = mat["ncounts"];
    IntegerVector diag_idx = mat["diag.idx"];
    IntegerVector diag_grp = mat["diag.grp"];
    NumericVector beta_cv = mat["beta_cv"];
    IntegerVector cv_grp = mat["cv.group"];
    //get patch nos and sorted values
    IntegerVector patchno = get_patch_numbers(nbins, mat, tol_val);
    //NumericVector patchvals = abs(get_patch_values(beta, patchno));
    NumericVector patchvals = abs(get_patch_values(beta_cv, patchno));
    std::sort(patchvals.begin(),patchvals.end());
    //since constraint is on, decay and signal must adjust so that
    //there is at least one zero signal value per diagonal idx
    NumericVector forbidden_vals;
    if (constrained) {
        if (positive) {
          forbidden_vals = get_minimum_diagonal_values(beta, diag_grp);
          lmin = std::max(lmin, std::max(std::abs(forbidden_vals(0)), std::abs(forbidden_vals(forbidden_vals.size()-1))));
        } else {
          NumericVector abeta=abs(beta);
          forbidden_vals = get_minimum_diagonal_values(abeta, diag_grp);
        }
    }
    //create functor
    /*obj_lambda1_BIC obj(lmin, tol_val, patchno, forbidden_vals,
                        beta, weight, phihat, ncounts);*/
    obj_lambda1_CV obj(lmin, tol_val, patchno, forbidden_vals,
                    beta_cv, weight, phihat, ncounts, cv_grp);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //get minimum authorized patch value
    double minpatch = max(abs(forbidden_vals));
    const double k=1;
    NumericVector best = optimize_CV_kSD(obj, 0, minpatch, max(abs(beta)) + 2*tol_val, tol_val, patchvals, k);
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["BIC.sd"]=best["BIC.sd"], _["dof"]=best["dof"]);
}
