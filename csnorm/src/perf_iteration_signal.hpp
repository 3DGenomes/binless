#ifndef PERF_ITERATION_SIGNAL_HPP
#define PERF_ITERATION_SIGNAL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           double lam1, double lam2,  double eCprime,
                           double alpha, double converge,
                           List outliers, NumericVector phi_i);

List wgfl_signal_cv(const DataFrame mat, int nbins, int ntrails,
                    double lam2, double alpha, double converge, NumericVector beta_i);

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int nbins,
                     double lam2,  double alpha, double tol_val,
                     List outliers, NumericVector beta_i, double lambda1_min, int refine_num,
                     bool constrained, bool fixed);

List wgfl_signal_BIC_fixed(const DataFrame cts, double dispersion, int nouter, int nbins,
                     double lam1, double lam2, double eCprime, double alpha, double tol_val,
                     List outliers, NumericVector beta_i);


#endif

