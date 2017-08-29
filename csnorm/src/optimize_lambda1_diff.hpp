#ifndef OPTIMIZE_LAMBDA1_DIFF_HPP
#define OPTIMIZE_LAMBDA1_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <utility> //pair

#include "optimize_lambda1.hpp" //obj_lambda1_base and compute_CV_diff
#include "Scores.hpp"
#include "ScoreComputer.hpp"
#include "DataLikelihoods.hpp"

//objective functor to find lambda1 assuming eCprime=0, using BIC
struct obj_lambda1_diff_BIC : private obj_lambda1_base,
                              private ScoreComputer<DifferenceLikelihood,BICScore> {
  obj_lambda1_diff_BIC(double minUB, double tol_val,
                  IntegerVector patchno, NumericVector forbidden_vals,
                  NumericVector value, NumericVector weight, NumericVector valuehat,
                  NumericVector weight_ref, NumericVector valuehat_ref,
                  NumericVector ncounts);
  
  double operator()(double x) const;
  
  NumericVector get(double val, std::string msg = "") const;
  
  double minUB_, tol_val_;
  NumericVector forbidden_vals_;
};

//objective functor to find lambda1 assuming eCprime=0, using CV
struct obj_lambda1_diff_CV : private obj_lambda1_base,
                             private ScoreComputer<DifferenceLikelihood,CVScore> {
    obj_lambda1_diff_CV(double minUB, double tol_val,
                IntegerVector patchno, NumericVector forbidden_vals,
                NumericVector value, NumericVector weight, NumericVector valuehat,
                NumericVector weight_ref, NumericVector valuehat_ref,
                NumericVector ncounts, IntegerVector cv_grp);

    double operator()(double x) const;

    NumericVector get(double val, std::string msg = "") const;

    double minUB_, tol_val_;
    NumericVector forbidden_vals_;
};

NumericVector cpp_optimize_lambda1_diff(const DataFrame mat, int nbins,
                                   double tol_val, double lambda1_min, int refine_num);


#endif

