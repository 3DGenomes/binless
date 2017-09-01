#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "optimize_lambda1.hpp"
#include "util.hpp"
#include "graph_helpers.hpp" //get_patch_numbers
#include "Traits.hpp"
#include "SparsityEstimator.hpp"



NumericVector cpp_optimize_lambda1(const DataFrame mat, int nbins,
                                   double tol_val, bool positive,
                                   double lambda1_min, int refine_num, double lambda2) {
    /*const bool constrained = true;
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
    SignalData data(beta, weight, phihat, ncounts, patchno); //TODO: beta or beta_cv ?
    //obj_lambda1<BIC> obj(lmin, tol_val, data, forbidden_vals, beta, ncounts);
    obj_lambda1<CV> obj(lmin, tol_val, data, forbidden_vals, beta_cv, cv_grp);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //get minimum authorized patch value
    double minpatch = max(abs(forbidden_vals));
    const double k=1;
    NumericVector best = optimize_CV_kSD(obj, 0, minpatch, max(abs(beta)) + 2*tol_val, tol_val, patchvals, k);
    //finalize
    //obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["BIC.sd"]=best["BIC.sd"], _["dof"]=best["dof"]);*/
    NumericVector beta = mat["beta"];
    NumericVector weight = mat["weight"];
    NumericVector phihat = mat["phihat"];
    NumericVector ncounts = mat["ncounts"];
    IntegerVector patchno = get_patch_numbers(nbins, mat, tol_val);
    NumericVector beta_cv = mat["beta_cv"];
    IntegerVector cv_grp = mat["cv.group"];
    SignalData data(beta, weight, phihat, ncounts, patchno);
    SparsityEstimator<Signal, CV, ZeroOffset, PositiveSign, ForbidDegeneracy> est(nbins, tol_val, data, lambda2, mat, beta_cv, cv_grp);
    return est.optimize();
}
