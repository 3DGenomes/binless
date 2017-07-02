#ifndef PERF_ITERATION_DIFF_HPP
#define PERF_ITERATION_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

DataFrame cts_to_diff_mat(const DataFrame cts, const DataFrame ref, int nbins,
                          double dispersion,
                          std::vector<double>& phi_ref, std::vector<double>& delta, int diag_rm);

std::vector<double> compute_phi_ref(const std::vector<double>& delta_r,
                                    const std::vector<double>& phihat,
                                    const std::vector<double>& phihat_var, const std::vector<double>& phihat_ref,
                                    const std::vector<double>& phihat_var_ref);

List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins,
                         int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                         double lam1, double lam2, double alpha, double inflate, int ninner,
                         double converge,
                         int diag_rm, NumericVector phi_ref_i, NumericVector beta_i);

List wgfl_diff_cv(const DataFrame mat, int nbins,
                    int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                    double lam2, double alpha, double inflate, int ninner, double converge, NumericVector beta_i);

List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                   double lam2,  double alpha, double inflate, int ninner, double tol_val,
                   int diag_rm, NumericVector phi_ref_i, NumericVector beta_i, double lambda1_min,
                   int refine_num, bool constrained);

List wgfl_diff_BIC_fixed(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                   double lam1, double lam2,  double alpha, double inflate, int ninner, double tol_val,
                   int diag_rm, NumericVector phi_ref_i, NumericVector beta_i);

#endif

