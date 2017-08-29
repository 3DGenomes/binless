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
                            beta, weight, phihat, lambda2, ncounts);*/
    obj_lambda1_eCprime<CVScore> obj(tol_val, constrained, patchno,
                               forbidden_vals, beta_cv, weight, phihat, lambda2, cv_grp);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    double minpatch = max(forbidden_vals);
    NumericVector best = optimize_CV(obj, patchvals(0) - 2*tol_val, minpatch, maxval + 2*tol_val, tol_val, patchvals);
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["BIC.sd"]=best["BIC.sd"], _["dof"]=best["dof"]);
}
