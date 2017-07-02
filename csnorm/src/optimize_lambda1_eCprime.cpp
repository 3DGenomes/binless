#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <string>

#include "optimize_lambda1_eCprime.hpp"
#include "util.hpp"
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include <boost/math/tools/minima.hpp> //brent_find_minima

obj_lambda1_eCprime_BIC::obj_lambda1_eCprime_BIC(double tol_val,
        bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
        NumericVector value, NumericVector weight, NumericVector valuehat,
        NumericVector ncounts, double lambda2) :
    tol_val_(tol_val), lsnc_(log(sum(ncounts))), lambda2_(lambda2), constrained_(constrained),
    patchno_(patchno), forbidden_vals_(forbidden_vals),
    value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_eCprime_BIC::operator()(double val) const {
    //return get(val, "opt")["BIC"];
    return get(val)["BIC"];
}

//take val as a starting point for optimization of UB and LB at constant dof
NumericVector obj_lambda1_eCprime_BIC::get(double val, std::string msg) const {
    //split data in two groups
    LogicalVector grp1 = value_ <= val;
    double b1,b2, xk, xkp1;
    double a,b,UB,LB;
    bool grp1only = is_true(all(grp1)), grp2only = is_true(all(!grp1));
    double xmin = min(value_), xmax = max(value_);
    
    //compute UB and LB
    if ((!grp1only) && (!grp2only)) {//non-degenerate case
      NumericVector w1 = weight_[grp1];
      NumericVector betahat1 = valuehat_[grp1];
      NumericVector w2 = weight_[!grp1];
      NumericVector betahat2 = valuehat_[!grp1];
      NumericVector x1 = value_[grp1];
      NumericVector x2 = value_[!grp1];
      b1 = sum(w1*betahat1)/sum(w1);
      b2 = sum(w2*(x2-betahat2))/sum(w2);
      xk = max(x1);
      xkp1 = min(x2);
      //compute unconstrained minimum for UB and LB
      a = b1 + b2;
      b = b1 - b2;
      //compute optimal UB and LB
      UB = std::max(std::min(xkp1, a), xk);
      LB = std::min(b, xmin);
      
    } else if (grp1only) {//grp2 is empty, UB > xmax
      NumericVector w1 = weight_[grp1];
      NumericVector betahat1 = valuehat_[grp1];
      b1 = sum(w1*betahat1)/sum(w1);
      b2 = 0;
      xk = xmax;
      xkp1 = std::numeric_limits<double>::infinity();
      //compute unconstrained minimum for UB and LB
      a = b1 + b2;
      b = b1 - b2;
      //compute optimal UB and LB
      if (b1 >= (xmin+xmax)/2.) {
        LB = xmin; //or any value < xmin
        UB = 2*b1 - LB;
      } else {
        UB = xmax; //or any value > xmax
        LB = 2*b1 - UB;
      }
      
    } else {//grp1 is empty, UB <= xmin
        NumericVector w2 = weight_[!grp1];
        NumericVector betahat2 = valuehat_[!grp1];
        NumericVector x2 = value_[!grp1];
        b1 = 0;
        b2 = sum(w2*(x2-betahat2))/sum(w2);
        xk = -std::numeric_limits<double>::infinity();
        xkp1 = xmin;
        //compute unconstrained minimum for UB and LB
        a = b1 + b2;
        b = b1 - b2;
        //compute optimal UB and LB
        UB = xmin;
        if (b2<-tol_val_) Rcout << " warning: b2= " << b2 << " is negative, setting UB=LB" << std::endl;
        LB = UB - 2*(std::max(b2,0.)+tol_val_);
    }
    
    //compute eCprime, lambda1
    double eCprime = (UB+LB)/2;
    double lambda1 = (UB-LB)/2;
    /*Rcout << "EVAL at val= " << val << " LB= " << LB << " UB= " << UB << " b1= " << b1 << " b2= " << b2
            << " xmin= " << xmin << " xk= " << xk << " xkp1= " << xkp1 << " a= " << a << " b= " << b << std::endl; */
    //check if solution is feasible
    if (constrained_) {
        if ( is_true(any( (forbidden_vals_>UB+tol_val_/2) | (forbidden_vals_<LB-tol_val_/2) )) ) {
            if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda2= " << lambda2_ << " lambda1= " << lambda1
              << " eCprime= " << eCprime << " BIC= Inf dof= NA" << std::endl;
            return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                         _["BIC"]=std::numeric_limits<double>::max(), _["dof"]=NumericVector::get_na(),
                                         _["UB"]=UB, _["LB"]=LB);
        }
    }
    //compute dof and BIC
    std::vector<double> value_r = as<std::vector<double> >(value_);
    NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
    IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
    const int dof = unique(selected).size();
    const double BIC = sum(weight_ * SQUARE(valuehat_ - (soft + eCprime))) +
                       lsnc_*dof;
    if (!msg.empty()) Rcout << " OBJ " << msg << " ok lambda2= " << lambda2_ << " lambda1= " << lambda1
                            << " eCprime= " << eCprime << " BIC= " << BIC  << " dof= " << dof 
                            << " UB= " << UB  << " LB= " << LB << std::endl;
    return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                 _["BIC"]=BIC, _["dof"]=dof, _["UB"]=UB, _["LB"]=LB);
}


obj_lambda1_eCprime_CV::obj_lambda1_eCprime_CV(double minval, double tol_val,
                                                 bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                                                 NumericVector value, NumericVector weight, NumericVector valuehat,
                                                 NumericVector ncounts, double lambda2) : minval_(minval),
                                                 tol_val_(tol_val), lsnc_(log(sum(ncounts))), lambda2_(lambda2), constrained_(constrained),
                                                 patchno_(patchno), forbidden_vals_(forbidden_vals),
                                                 value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_eCprime_CV::operator()(double val) const {
  //return get(val, "opt")["CV"];
  return get(val)["CV"];
}

//take val as a starting point for optimization of UB and LB at constant dof
NumericVector obj_lambda1_eCprime_CV::get(double val, std::string msg) const {
  //split data in two groups
  LogicalVector grp1 = value_ <= val;
  double b1,b2, xk, xkp1;
  double a,b,UB,LB;
  bool grp1only = is_true(all(grp1)), grp2only = is_true(all(!grp1));
  double xmin = min(value_), xmax = max(value_); //here, we assume xmin >= minval_
  
  //compute UB and LB
  if ((!grp1only) && (!grp2only)) {//non-degenerate case
    NumericVector w1 = weight_[grp1];
    NumericVector betahat1 = valuehat_[grp1];
    NumericVector w2 = weight_[!grp1];
    NumericVector betahat2 = valuehat_[!grp1];
    NumericVector x1 = value_[grp1];
    NumericVector x2 = value_[!grp1];
    b1 = sum(w1*betahat1)/sum(w1);
    b2 = sum(w2*(x2-betahat2))/sum(w2);
    xk = max(x1);
    xkp1 = min(x2);
    //compute unconstrained minimum for UB and LB
    a = b1 + b2;
    b = b1 - b2;
    //compute optimal UB and LB
    UB = std::max(std::min(xkp1, a), xk);
    LB = std::min(b, minval_);
    
  } else if (grp1only) {//grp2 is empty, UB > xmax
    NumericVector w1 = weight_;
    NumericVector betahat1 = valuehat_;
    b1 = sum(w1*betahat1)/sum(w1);
    b2 = 0;
    xk = xmax;
    xkp1 = std::numeric_limits<double>::infinity();
    //compute unconstrained minimum for UB and LB
    a = b1 + b2;
    b = b1 - b2;
    //compute optimal UB and LB
    if (b1 >= (minval_+xmax)/2.) {
      LB = minval_; //or any value < minval_
      UB = 2*b1 - LB;
    } else {
      UB = xmax; //or any value > xmax
      LB = 2*b1 - UB;
    }
    
  } else {//grp1 is empty, UB <= xmin
    NumericVector w2 = weight_;
    NumericVector betahat2 = valuehat_;
    NumericVector x2 = value_;
    b1 = 0;
    b2 = sum(w2*(x2-betahat2))/sum(w2);
    xk = -std::numeric_limits<double>::infinity();
    xkp1 = xmin;
    //compute unconstrained minimum for UB and LB
    a = b1 + b2;
    b = b1 - b2;
    //compute optimal UB and LB
    if (b2 >= (xmin-minval_)/2.) {
      UB = xmin;
      LB = UB - 2*b2;
    } else {
      LB = minval_;
      if (b2 < -tol_val_) Rcout << " warning: b2= " << b2 << " is negative, setting UB=LB" << std::endl;
      UB = LB + 2*(std::max(b2,0.)+tol_val_);
    }
  }
  
  //compute eCprime, lambda1
  double eCprime = (UB+LB)/2;
  double lambda1 = (UB-LB)/2;
  /*Rcout << "EVAL at val= " << val << " LB= " << LB << " UB= " << UB << " b1= " << b1 << " b2= " << b2
          << " xmin= " << xmin << " xk= " << xk << " xkp1= " << xkp1 << " a= " << a << " b= " << b << std::endl; */
  //check if solution is feasible
  if (constrained_) {
    if ( is_true(any( (forbidden_vals_>UB+tol_val_/2) | (forbidden_vals_<LB-tol_val_/2) )) ) {
      if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda2= " << lambda2_ << " lambda1= " << lambda1
                              << " eCprime= " << eCprime << " CV= Inf dof= NA" << std::endl;
      return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                   _["BIC"]=std::numeric_limits<double>::max(), _["dof"]=NumericVector::get_na(),
                                   _["UB"]=UB, _["LB"]=LB);
    }
  }
  //compute dof and CV
  std::vector<double> value_r = as<std::vector<double> >(value_);
  NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
  IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
  const int dof = unique(selected).size();
  const double CV = mean(weight_ * SQUARE(valuehat_ - (soft + eCprime)));
  if (!msg.empty()) Rcout << " OBJ " << msg << " ok lambda2= " << lambda2_ << " lambda1= " << lambda1
                          << " eCprime= " << eCprime << " CV= " << CV  << " dof= " << dof 
                          << " UB= " << UB  << " LB= " << LB << std::endl;
  return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                               _["BIC"]=CV, _["dof"]=dof, _["UB"]=UB, _["LB"]=LB);
}


NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num, double lambda2) {
    //extract vectors
    double lmin = std::max(lambda1_min,tol_val/2);
    std::clock_t c_start = std::clock();
    NumericVector weight = mat["weight"];
    NumericVector phihat = mat["phihat"];
    NumericVector beta = mat["beta"];
    NumericVector ncounts = mat["ncounts"];
    IntegerVector diag_idx = mat["diag.idx"];
    NumericVector beta_cv = mat["beta_cv"];
    //get patch nos and sorted values
    List cl = boost_build_patch_graph_components(nbins, mat, tol_val);
    IntegerVector patchno = cl["membership"];
    //NumericVector patchvals = get_patch_values(beta, patchno);
    NumericVector patchvals = get_patch_values(beta_cv, patchno);
    double minval = min(beta);
    //if constraint is on, decay and signal must adjust so that
    //there is at least one zero signal value per diagonal idx
    NumericVector forbidden_vals;
    if (constrained) {
      forbidden_vals = get_minimum_diagonal_values(beta, diag_idx);
      lmin = std::max(lmin, (max(forbidden_vals)-minval)/2);
    }
    //create functor
    /*obj_lambda1_eCprime_BIC obj(tol_val, constrained, patchno,
                            forbidden_vals,
                            beta, weight, phihat, ncounts, lambda2);*/
    obj_lambda1_eCprime_CV obj(minval, tol_val, constrained, patchno,
                               forbidden_vals,
                               beta_cv, weight, phihat, ncounts, lambda2);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //loop over patch values
    std::clock_t c_in1 = std::clock();
    //NumericVector best = obj.get(patchvals(0) - 2*tol_val, "opt"); //query is < xmin
    NumericVector best = obj.get(patchvals(0) - 2*tol_val); //query is < xmin
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
