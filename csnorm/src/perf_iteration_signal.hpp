#ifndef PERF_ITERATION_SIGNAL_HPP
#define PERF_ITERATION_SIGNAL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins, List GFLState,
                           double lam2, double converge,
                           const List outliers, NumericVector phi_i);

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int nbins, List GFLState,
                     double lam2, double tol_val,
                     List outliers, NumericVector beta_i, bool constrained, bool fixed);


#endif

