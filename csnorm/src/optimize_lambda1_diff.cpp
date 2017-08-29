#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "optimize_lambda1_diff.hpp"
#include "util.hpp"
#include "graph_helpers.hpp" //get_patch_numbers

obj_lambda1_diff_BIC::obj_lambda1_diff_BIC(double minUB, double tol_val,
                                 IntegerVector patchno, NumericVector forbidden_vals,
                                 NumericVector value, NumericVector weight, NumericVector valuehat,
                                 NumericVector weight_ref, NumericVector valuehat_ref,
                                 NumericVector ncounts) :
  obj_lambda1_base(value, weight, valuehat, minUB),
  ScoreComputer(tol_val, value, weight, valuehat, weight_ref, valuehat_ref, patchno, ncounts),
  minUB_(minUB), tol_val_(tol_val), lsnc_(log(sum(ncounts))), forbidden_vals_(forbidden_vals),
  patchno_(patchno), value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_diff_BIC::operator()(double x) const {
  //return get(std::pow(10,x), "opt")["BIC"];
  return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_diff_BIC::get(double val, std::string msg) const {
  
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
  std::vector<double> value_r = as<std::vector<double> >(value_);
  NumericVector soft = wrap(soft_threshold(value_r, 0, lambda1));
  IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
  const int dof = unique(selected).size();
  const double BIC = sum(weight_ * SQUARE(valuehat_ - soft)) + lsnc_*dof;
  if (!msg.empty()) Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
                          << " BIC= " << BIC  << " dof= " << dof
                          << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
  return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1, _["BIC"]=BIC,
                               _["dof"]=dof, _["UB"]=UB, _["LB"]=-UB);
}

obj_lambda1_diff_CV::obj_lambda1_diff_CV(double minUB, double tol_val,
                         IntegerVector patchno, NumericVector forbidden_vals,
                         NumericVector value, NumericVector weight, NumericVector valuehat,
                         NumericVector weight_ref, NumericVector valuehat_ref,
                         NumericVector ncounts, IntegerVector cv_grp) :
    obj_lambda1_base(value, weight, valuehat, minUB),
    ScoreComputer(tol_val, value, weight, valuehat, weight_ref, valuehat_ref, patchno, cv_grp),
    minUB_(minUB), tol_val_(tol_val), lsnc_(log(sum(ncounts))), forbidden_vals_(forbidden_vals),
    patchno_(patchno), value_(value), weight_(weight), valuehat_(valuehat), weight_ref_(weight_ref),
    valuehat_ref_(valuehat_ref), cv_grp_(cv_grp) {}

double obj_lambda1_diff_CV::operator()(double x) const {
    //return get(std::pow(10,x), "opt")["BIC"];
    return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_diff_CV::get(double val, std::string msg) const {
  
    
    double UB = optimize_bounds(val);
    
    
    //check if forbidden
  double lambda1=UB;
  if ( is_true(any(abs(forbidden_vals_)>lambda1+tol_val_/2)) || (UB < minUB_) ) {
      if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1
                              << " eCprime= 0 CV= Inf dof= NA"
                              << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
      return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                   _["BIC"]=std::numeric_limits<double>::max(), _["BIC.sd"]=0,
                                   _["dof"]=NumericVector::get_na(),
                                   _["UB"]=lambda1, _["LB"]=-lambda1);
  }
    return evaluate(-UB, UB);
}

NumericVector cpp_optimize_lambda1_diff(const DataFrame mat, int nbins,
                                   double tol_val,
                                   double lambda1_min, int refine_num) {
    const bool constrained = true;
    //extract vectors
    double lmin = std::max(lambda1_min,tol_val/2);
    NumericVector weight = 1/as<NumericVector>(mat["phihat.var"]);
    NumericVector weight_ref = 1/as<NumericVector>(mat["phihat.var.ref"]);
    NumericVector phihat = mat["phihat"];
    NumericVector phihat_ref = mat["phihat.ref"];
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
      NumericVector abeta=abs(beta);
      forbidden_vals = get_minimum_diagonal_values(abeta, diag_grp);
    }
    //create functor
    /*obj_lambda1_diff_BIC obj(lmin, tol_val, patchno, forbidden_vals,
                        beta, weight, phihat, ncounts);*/
    obj_lambda1_diff_CV obj(lmin, tol_val, patchno, forbidden_vals,
                    beta_cv, weight, phihat, weight_ref, phihat_ref, ncounts, cv_grp);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //get minimum authorized patch value
    double minpatch = max(abs(forbidden_vals));
    const double k=1;
    NumericVector best = optimize_CV_kSD(obj, 0, minpatch, 2*max(abs(beta)) + 2*tol_val, tol_val, patchvals, k);
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["BIC.sd"]=best["BIC.sd"], _["dof"]=best["dof"]);
}
