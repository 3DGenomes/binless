#ifndef PERF_ITERATION_H
#define PERF_ITERATION_H

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

DataFrame cts_to_signal_mat(const DataFrame cts, int nbins, double dispersion, std::vector<double>& phi, int diag_rm);

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                           double lam,  double alpha, double inflate, int ninner, double converge,
                           int diag_rm, NumericVector z_i, NumericVector u_i, NumericVector phi_i);

List wgfl_signal_perf(const DataFrame cts, double dispersion, int nouter, int nbins,
                      int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                      double lam,  double alpha, double inflate, int ninner, double converge, int diag_rm);

DataFrame cts_to_diff_mat(const DataFrame cts, const DataFrame ref, int nbins, double dispersion,
                          std::vector<double>& phi_ref, std::vector<double>& delta, int diag_rm);

List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref, double dispersion, int nouter, int nbins,
                         int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                         double lam,  double alpha, double inflate, int ninner, double converge,
                         int diag_rm, NumericVector z_i, NumericVector u_i, NumericVector phi_ref_i,
                         NumericVector delta_i);

List wgfl_diff_perf(const DataFrame cts, const DataFrame ref, double dispersion, int nouter, int nbins,
                    int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                    double lam,  double alpha, double inflate, int ninner, double converge, int diag_rm);

#endif
