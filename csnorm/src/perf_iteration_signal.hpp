#ifndef PERF_ITERATION_SIGNAL_HPP
#define PERF_ITERATION_SIGNAL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int nbins, List GFLState,
                     double lam2, double tol_val,
                     List metadata, NumericVector beta_i, bool constrained, bool fixed);


#endif

