#ifndef OPTIMIZE_LAMBDA1_DIFF_HPP
#define OPTIMIZE_LAMBDA1_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <utility> //pair

#include "optimize_lambda1.hpp" //obj_lambda1_base

struct compute_CV_diff {
    compute_CV_diff(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
                    const NumericVector& weight_ref, const NumericVector& valuehat_ref,
                    const IntegerVector& patchno, const IntegerVector& cv_grp) :
    tol_val_(tol_val), value_(value), weight_(weight), valuehat_(valuehat), weight_ref_(weight_ref),
    valuehat_ref_(valuehat_ref), patchno_(patchno), cv_grp_(cv_grp) {}
    NumericVector evaluate(double LB, double UB) const;
private:
    double tol_val_;
    NumericVector value_, weight_, valuehat_, weight_ref_, valuehat_ref_;
    IntegerVector patchno_, cv_grp_;
};

//objective functor to find lambda1 assuming eCprime=0, using BIC
struct obj_lambda1_diff_BIC : private obj_lambda1_base {
  obj_lambda1_diff_BIC(double minUB, double tol_val,
                  IntegerVector patchno, NumericVector forbidden_vals,
                  NumericVector value, NumericVector weight, NumericVector valuehat,
                  NumericVector ncounts);
  
  double operator()(double x) const;
  
  NumericVector get(double val, std::string msg = "") const;
  
  double minUB_, tol_val_, lsnc_;
  NumericVector forbidden_vals_;
  IntegerVector patchno_;
  NumericVector value_, weight_, valuehat_;
};

//objective functor to find lambda1 assuming eCprime=0, using CV
struct obj_lambda1_diff_CV : private obj_lambda1_base, private compute_CV_diff {
    obj_lambda1_diff_CV(double minUB, double tol_val,
                IntegerVector patchno, NumericVector forbidden_vals,
                NumericVector value, NumericVector weight, NumericVector valuehat,
                NumericVector weight_ref, NumericVector valuehat_ref,
                NumericVector ncounts, IntegerVector cv_grp);

    double operator()(double x) const;

    NumericVector get(double val, std::string msg = "") const;

    double minUB_, tol_val_, lsnc_;
    NumericVector forbidden_vals_;
    IntegerVector patchno_;
    NumericVector value_, weight_, valuehat_, weight_ref_, valuehat_ref_;
    IntegerVector cv_grp_;
};

NumericVector cpp_optimize_lambda1_diff(const DataFrame mat, int nbins,
                                   double tol_val, double lambda1_min, int refine_num);


#endif

