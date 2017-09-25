#ifndef PERF_ITERATION_DIFF_HPP
#define PERF_ITERATION_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins, List GFLState,
                         double lam2, double converge,
                         List outliers, NumericVector phi_ref_i, NumericVector beta_i);

List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins, List GFLState,
                   double lam2,  double tol_val,
                   List outliers, NumericVector phi_ref_i, NumericVector beta_i, double lambda1_min,
                   int refine_num, bool constrained);

#endif

