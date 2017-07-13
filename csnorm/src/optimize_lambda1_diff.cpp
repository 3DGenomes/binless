#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <string>

#include "optimize_lambda1_diff.hpp"
#include "util.hpp"
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include <boost/math/tools/minima.hpp> //brent_find_minima

obj_lambda1_diff_BIC::obj_lambda1_diff_BIC(double minUB, double tol_val,
                                 IntegerVector patchno, NumericVector forbidden_vals,
                                 NumericVector value, NumericVector weight, NumericVector valuehat,
                                 NumericVector ncounts) :
  minUB_(minUB), minabsval_(min(abs(value))), maxabsval_(max(abs(value))),
  tol_val_(tol_val), lsnc_(log(sum(ncounts))),
  patchno_(patchno), forbidden_vals_(forbidden_vals),
  absval_(abs(value)), value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_diff_BIC::operator()(double x) const {
  //return get(std::pow(10,x), "opt")["BIC"];
  return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_diff_BIC::get(double val, std::string msg) const {
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
  if (minUB_ <= xkp1) UB=std::max(minUB_,UB);
  //Rcout << "  UB= " << UB << " minUB= " << minUB_ << std::endl;
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
                         NumericVector ncounts) :
    minUB_(minUB), minabsval_(min(abs(value))), maxabsval_(max(abs(value))),
    tol_val_(tol_val), lsnc_(log(sum(ncounts))),
    patchno_(patchno), forbidden_vals_(forbidden_vals),
    absval_(abs(value)), value_(value), weight_(weight), valuehat_(valuehat),
    weight_ref_(weight_ref), valuehat_ref_(valuehat_ref) {}

double obj_lambda1_diff_CV::operator()(double x) const {
    //return get(std::pow(10,x), "opt")["BIC"];
    return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_diff_CV::get(double val, std::string msg) const {
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
  //check if forbidden
  double lambda1=UB;
  if ( is_true(any(abs(forbidden_vals_)>lambda1+tol_val_/2)) || (UB < minUB_) ) {
      if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1
                              << " eCprime= 0 CV= Inf dof= NA"
                              << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
      return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                   _["BIC"]=std::numeric_limits<double>::max(), _["dof"]=NumericVector::get_na(),
                                   _["UB"]=lambda1, _["LB"]=-lambda1);
  }
    std::vector<double> value_r = as<std::vector<double> >(value_);
    std::vector<double> soft_r = soft_threshold(value_r, 0, lambda1);
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
    const double CV = sum(weight_ref_ * SQUARE(valuehat_ref_ - phi_ref) + weight_*SQUARE(valuehat_-(phi_ref+soft)));
    if (!msg.empty()) Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
                            << " CV= " << CV  << " dof= " << dof
                            << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
    return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1, _["BIC"]=CV,
                                 _["dof"]=dof, _["UB"]=UB, _["LB"]=-UB);
}

NumericVector cpp_optimize_lambda1_diff(const DataFrame mat, int nbins,
                                   double tol_val,
                                   double lambda1_min, int refine_num) {
    const bool constrained = true;
    //extract vectors
    double lmin = std::max(lambda1_min,tol_val/2);
    std::clock_t c_start = std::clock();
    NumericVector weight = 1/as<NumericVector>(mat["phihat.var"]);
    NumericVector weight_ref = 1/as<NumericVector>(mat["phihat.var.ref"]);
    NumericVector phihat = mat["phihat"];
    NumericVector phihat_ref = mat["phihat.ref"];
    NumericVector beta = mat["beta"];
    NumericVector ncounts = mat["ncounts"];
    IntegerVector diag_idx = mat["diag.idx"];
    NumericVector beta_cv = mat["beta_cv"];
    //get patch nos and sorted values
    List cl = boost_build_patch_graph_components(nbins, mat, tol_val);
    IntegerVector patchno = cl["membership"];
    //NumericVector patchvals = abs(get_patch_values(beta, patchno));
    NumericVector patchvals = abs(get_patch_values(beta_cv, patchno));
    std::sort(patchvals.begin(),patchvals.end());
    //since constraint is on, decay and signal must adjust so that
    //there is at least one zero signal value per diagonal idx
    NumericVector forbidden_vals;
    if (constrained) {
      NumericVector abeta=abs(beta);
      forbidden_vals = get_minimum_diagonal_values(abeta, diag_idx);
    }
    //create functor
    /*obj_lambda1_diff_BIC obj(lmin, tol_val, patchno, forbidden_vals,
                        beta, weight, phihat, ncounts);*/
    obj_lambda1_diff_CV obj(lmin, tol_val, patchno, forbidden_vals,
                    beta_cv, weight, phihat, weight_ref, phihat_ref, ncounts);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //get minimum authorized patch value
    double minpatch = max(abs(forbidden_vals));
    /*Rcout << "minpatch= " << minpatch << " npatches= " << patchvals.size() << " npatches.ok= "
            << as<NumericVector>(patchvals[patchvals>=minpatch]).size() << std::endl;*/
    //loop over patch values
    std::clock_t c_in1 = std::clock();
    //NumericVector best = obj.get(patchvals(0) - 2*tol_val, "opt"); //query is < xmin
    NumericVector best = obj.get(0); //query is < xmin
    for (int i=0; i<patchvals.size(); ++i) {
      if (patchvals(i) <= minpatch) continue;
      //NumericVector val = obj.get(patchvals(i) + 2*tol_val, "opt");
      NumericVector val = obj.get(patchvals(i) + 2*tol_val);
      if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
    }
    //NumericVector val = obj.get(max(abs(beta)) + 2*tol_val, "opt");
    NumericVector val = obj.get(max(abs(beta)) + 2*tol_val);
    if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
    std::clock_t c_in2 = std::clock();
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["dof"]=best["dof"],
                                 _["c_init"]=c_in1-c_start, _["c_brent"]=0, _["c_refine"]=c_in2-c_in1);
}