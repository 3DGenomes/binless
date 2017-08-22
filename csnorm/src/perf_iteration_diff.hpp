#ifndef PERF_ITERATION_DIFF_HPP
#define PERF_ITERATION_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

DataFrame cts_to_diff_mat(const DataFrame cts, const DataFrame ref, int nbins,
                          double dispersion,
                          std::vector<double>& phi_ref, std::vector<double>& delta, List outliers);

List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins,
                         double lam1, double lam2, double alpha,
                         double converge,
                         List outliers, NumericVector phi_ref_i, NumericVector beta_i);

List wgfl_diff_cv(const DataFrame mat, int nbins,
                    double lam2, double alpha, double converge, NumericVector beta_i);

List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   double lam2,  double alpha, double tol_val,
                   List outliers, NumericVector phi_ref_i, NumericVector beta_i, double lambda1_min,
                   int refine_num, bool constrained);

List wgfl_diff_BIC_fixed(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   double lam1, double lam2,  double alpha, double tol_val,
                   List outliers, NumericVector phi_ref_i, NumericVector beta_i);

#endif

