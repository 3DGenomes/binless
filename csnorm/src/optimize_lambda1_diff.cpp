#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "optimize_lambda1_diff.hpp"
#include "util.hpp"
#include "graph_helpers.hpp" //get_patch_numbers
#include "Scores.hpp"

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
    /*obj_lambda1_diff<BICScore> obj(lmin, tol_val, patchno, forbidden_vals,
                        beta, weight, phihat, ncounts, ncounts);*/
    obj_lambda1_diff<CVScore> obj(lmin, tol_val, patchno, forbidden_vals,
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
