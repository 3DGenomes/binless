#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <string>

#include "optimize_lambda1.hpp"
#include "util.hpp"
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include <boost/math/tools/minima.hpp> //brent_find_minima

obj_lambda1_BIC::obj_lambda1_BIC(double minUB, double tol_val,
                         IntegerVector patchno, NumericVector forbidden_vals,
                         NumericVector value, NumericVector weight, NumericVector valuehat,
                         NumericVector ncounts) :
    minUB_(minUB), minabsval_(min(abs(value_))), maxabsval_(max(abs(value_))),
    tol_val_(tol_val), lsnc_(log(sum(ncounts))),
    patchno_(patchno), forbidden_vals_(forbidden_vals),
    absval_(abs(value)), value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_BIC::operator()(double x) const {
    //return get(std::pow(10,x), "opt")["BIC"];
    return get(std::pow(10,x))["BIC"];
}

NumericVector obj_lambda1_BIC::get(double val, std::string msg) const {
  //split data in two groups and determine constraint values
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
  //determine unconstrained UB
  double UB;
  bool grp1empty = is_true(all(!grp1)), grp2empty = is_true(all(!grp2));
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

NumericVector cpp_optimize_lambda1(const DataFrame mat, int nbins,
                                   double tol_val, bool positive,
                                   double lambda1_min, int refine_num) {
    const bool constrained = true;
    //extract vectors
    double lmin = std::max(lambda1_min,tol_val/2);
    std::clock_t c_start = std::clock();
    NumericVector weight = mat["weight"];
    NumericVector phihat = mat["phihat"];
    NumericVector beta = mat["beta"];
    NumericVector ncounts = mat["ncounts"];
    IntegerVector diag_idx = mat["diag.idx"];
    //get patch nos and sorted values
    List cl = boost_build_patch_graph_components(nbins, mat, tol_val);
    IntegerVector patchno = cl["membership"];
    NumericVector patchvals = abs(get_patch_values(beta, patchno));
    std::sort(patchvals.begin(),patchvals.end());
    //since constraint is on, decay and signal must adjust so that
    //there is at least one zero signal value per diagonal idx
    NumericVector forbidden_vals;
    if (constrained) {
        if (positive) {
            forbidden_vals = get_minimum_diagonal_values(beta, diag_idx);
            lmin = std::max(lmin, std::max(std::abs(forbidden_vals(0)), std::abs(forbidden_vals(forbidden_vals.size()-1))));
        } else {
            forbidden_vals = get_constant_diagonal_values(beta, diag_idx, tol_val);
        }
    }
    //create functor
    obj_lambda1_BIC obj(lmin, tol_val, patchno, forbidden_vals,
                    beta, weight, phihat, ncounts);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //loop over patch values
    std::clock_t c_in1 = std::clock();
    //NumericVector best = obj.get(patchvals(0) - 2*tol_val, "opt"); //query is < xmin
    NumericVector best = obj.get(0); //query is < xmin
    for (int i=0; i<patchvals.size(); ++i) {
      if (patchvals(i) <= max(forbidden_vals)) continue;
      //NumericVector val = obj.get(patchvals(i) + 2*tol_val, "opt");
      NumericVector val = obj.get(patchvals(i) + 2*tol_val);
      if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
    }
    std::clock_t c_in2 = std::clock();
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["dof"]=best["dof"],
                                 _["c_init"]=c_in1-c_start, _["c_brent"]=0, _["c_refine"]=c_in2-c_in1);
}
