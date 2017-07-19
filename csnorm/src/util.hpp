#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#define SQUARE(x) ((x)*(x))

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1);

NumericVector get_patch_values(NumericVector value, IntegerVector patchno);

NumericVector get_minimum_diagonal_values(NumericVector value,
        IntegerVector diag_idx);

NumericVector get_constant_diagonal_values(NumericVector value,
        IntegerVector diag_idx, double tol_val);

std::vector<double> compute_phi_ref(const std::vector<double>& delta_r,
                                    const std::vector<double>& phihat,
                                    const std::vector<double>& phihat_var, const std::vector<double>& phihat_ref,
                                    const std::vector<double>& phihat_var_ref);

//find minimum of obj in patchvals between startval and maxval, also evaluating minval
template<class Objective>
NumericVector optimize_CV(Objective obj, double minval, double startval, double maxval, double tol_val, NumericVector patchvals) {
  /*Rcout << "minval= " << minval << " startval= " << startval << " npatches= " << patchvals.size() << " npatches.ok= "
          << as<NumericVector>(patchvals[patchvals>=minval]).size() << std::endl;*/
  //loop over patch values
  //NumericVector best = obj.get(minval, "opt"); //query is < xmin
  NumericVector best = obj.get(minval); //query is < xmin
  for (int i=0; i<patchvals.size(); ++i) {
    if (patchvals(i) <= startval) continue;
    //NumericVector val = obj.get(patchvals(i) + 2*tol_val, "opt");
    NumericVector val = obj.get(patchvals(i) + 2*tol_val);
    if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
  }
  //NumericVector val = obj.get(maxval, "opt");
  NumericVector val = obj.get(maxval);
  if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
  return(best);
}

#endif

