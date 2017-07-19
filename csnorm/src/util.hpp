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

//evaluate obj at minval, maxval and all patchvals+2*tol_val > startval+2*tol_val
template<class Objective>
std::vector<NumericVector> optimize_CV_core(Objective obj, double minval, double startval, double maxval, double tol_val, NumericVector patchvals) {
  std::vector<NumericVector> ret;
 /*Rcout << "minval= " << minval << " startval= " << startval << " npatches= " << patchvals.size() << " npatches.ok= "
          << as<NumericVector>(patchvals[patchvals>=minval]).size() << std::endl;*/
  //loop over patch values
  //ret.push_back(obj.get(minval, "opt")); //query is < xmin
  ret.push_back(obj.get(minval)); //query is < xmin
  for (int i=0; i<patchvals.size(); ++i) {
    if (patchvals(i) <= startval) continue;
    //ret.push_back(obj.get(patchvals(i) + 2*tol_val, "opt"));
    ret.push_back(obj.get(patchvals(i) + 2*tol_val));
  }
  //ret.push_back(obj.get(maxval, "opt"));
  ret.push_back(obj.get(maxval));
  return(ret);
}

//get minimum of obj in patchvals between startval and maxval, also evaluating minval
template<class Objective>
NumericVector optimize_CV(Objective obj, double minval, double startval, double maxval, double tol_val, NumericVector patchvals) {
  std::vector<NumericVector> ret = optimize_CV_core(obj,minval,startval,maxval,tol_val,patchvals);
  NumericVector best = ret[0];
  for (unsigned i=1; i<ret.size(); ++i)
    if (as<double>(ret[i]["BIC"]) < as<double>(best["BIC"])) best=ret[i];
  return(best);
}

//get minimum+kSD of obj in patchvals between startval and maxval, also evaluating minval
template<class Objective>
NumericVector optimize_CV_kSD(Objective obj, double minval, double startval, double maxval, double tol_val, NumericVector patchvals, double k) {
  //evaluate everywhere
  std::vector<NumericVector> ret = optimize_CV_core(obj,minval,startval,maxval,tol_val,patchvals);
  //find global minimum
  NumericVector best = ret[0];
  unsigned best_idx = 0;
  for (unsigned i=1; i<ret.size(); ++i) {
    if (as<double>(ret[i]["BIC"]) < as<double>(best["BIC"])) {
      best=ret[i];
      best_idx=i;
    }
  }
  //find kSD solution
  double thresh = best["BIC"]+k*best["BIC.sd"];
  /*Rcout << " found optimum at " << as<double>(best["BIC"]) << " adding " << (k*as<double>(best["BIC.sd"]))
          << " thresh= " << thresh << " ret.size()= " << ret.size() << " best_idx= " << best_idx << " " << std::endl;*/
  unsigned i_opt = best_idx;
  for (i_opt = best_idx+1; i_opt < ret.size(); ++i_opt) {
      if (as<double>(ret[i_opt]["BIC"]) > thresh) break;
  }
  if (i_opt >= ret.size()) i_opt = ret.size()-1;
  return(ret[i_opt]);
}

#endif

