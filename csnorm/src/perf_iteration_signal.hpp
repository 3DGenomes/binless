#ifndef PERF_ITERATION_SIGNAL_HPP
#define PERF_ITERATION_SIGNAL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

DataFrame cts_to_signal_mat(const DataFrame cts, int nbins, double dispersion, std::vector<double>& phi,
                            double eCprime, int diag_rm);

std::vector<double> soft_threshold(const std::vector<double>& beta, double eCprime, double lam1);

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                           double lam1, double lam2,  double eCprime,
                           double alpha, double inflate, int ninner, double converge,
                           int diag_rm, NumericVector z_i, NumericVector u_i, NumericVector phi_i);

List wgfl_signal_perf(const DataFrame cts, double dispersion, int nouter, int nbins,
                      int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                      double lam1, double lam2, double eCprime,
                      double alpha, double inflate, int ninner, double converge, int diag_rm);
#endif

