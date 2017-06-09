#ifndef PERF_ITERATION_SIGNAL_HPP
#define PERF_ITERATION_SIGNAL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

DataFrame cts_to_signal_mat(const DataFrame cts, int nbins, double dispersion, std::vector<double>& phi,
                            double eCprime, int diag_rm);

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                           double lam1, double lam2,  double eCprime,
                           double alpha, double inflate, int ninner, double converge,
                           int diag_rm, NumericVector phi_i);

List wgfl_signal_perf_opt_lambda1_eCprime(const DataFrame cts, double dispersion, int nouter, int opt_every, int nbins,
                                          int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                                          double lam2, double alpha, double inflate, int ninner, double converge,
                                          int diag_rm,  double lam1_init, double eCprime_init, NumericVector beta_i,
                                          double lambda1_min, int refine_num);

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int opt_every, int nbins,
                     int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                     double lam2,  double alpha, double inflate, int ninner, double tol_val,
                     int diag_rm, NumericVector beta_i, double lambda1_min, int refine_num, bool constrained, bool fixed);


#endif

