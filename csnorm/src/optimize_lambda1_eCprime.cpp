#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include <utility> //pair
#include <tuple> //tie

#include "optimize_lambda1_eCprime.hpp"
#include "util.hpp"
#include "graph_helpers.hpp" //get_patch_numbers
#include "Scores.hpp"

    
obj_lambda1_eCprime_base::bounds_t obj_lambda1_eCprime_base::optimize_bounds(double val) const {
    //split data in two groups
    LogicalVector grp1 = value_ <= val;
    NumericVector w1 = weight_[grp1];
    NumericVector w2 = weight_[!grp1];
    bool grp1only = w2.size() == 0 || sum(w2) == 0;
    bool grp2only = w1.size() == 0 || sum(w1) == 0;
    double b1,b2, xk, xkp1;
    double a,b,UB,LB;
    double xmin = min(value_), xmax = max(value_); //here, we assume xmin >= minval_
  
  //compute UB and LB
  if ((!grp1only) && (!grp2only)) {//non-degenerate case
      NumericVector betahat1 = valuehat_[grp1];
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
      if (!grp2only) Rcout << "This should not have happened!\n";
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
          UB = LB + 2*b2;
      }
  }
  return bounds_t{LB, UB};
}

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num, double lambda2) {
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
    //NumericVector patchvals = get_patch_values(beta, patchno);
    NumericVector patchvals = get_patch_values(beta_cv, patchno);
    double minval = min(beta);
    double maxval = max(beta);
    //if constraint is on, decay and signal must adjust so that
    //there is at least one zero signal value per diagonal idx
    NumericVector forbidden_vals;
    if (constrained) {
      forbidden_vals = get_minimum_diagonal_values(beta, diag_grp);
      lmin = std::max(lmin, (max(forbidden_vals)-minval)/2);
    }
    //create functor
    /*obj_lambda1_eCprime<BICScore> obj(tol_val, constrained, patchno,
                            forbidden_vals,
                            beta, weight, phihat, ncounts, lambda2, ncounts);*/
    obj_lambda1_eCprime<CVScore> obj(tol_val, constrained, patchno,
                               forbidden_vals, beta_cv, weight, phihat, ncounts, lambda2, cv_grp);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    double minpatch = max(forbidden_vals);
    NumericVector best = optimize_CV(obj, patchvals(0) - 2*tol_val, minpatch, maxval + 2*tol_val, tol_val, patchvals);
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["BIC.sd"]=best["BIC.sd"], _["dof"]=best["dof"]);
}
