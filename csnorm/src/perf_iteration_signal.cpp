#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>

#include "perf_iteration_signal.hpp"
#include "gfl_c.h"
#include "graph_fl.h"

DataFrame cts_to_signal_mat(const DataFrame cts, int nbins, double dispersion, std::vector<double>& phi,
                            double eCprime, int diag_rm)
{
  //inputs
  int N = cts.nrows();
  std::vector<int> cts_bin1 = as<std::vector<int> >(cts["bin1"]);
  std::vector<int> cts_bin2 = as<std::vector<int> >(cts["bin2"]);
  std::vector<double> count = as<std::vector<double> >(cts["count"]);
  std::vector<double> lmu_nosig = as<std::vector<double> >(cts["lmu.nosig"]);
  std::vector<double> weight = as<std::vector<double> >(cts["weight"]);
  
  //outputs
  int nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<double> phihat(nbetas, 0); //vectorized form
  std::vector<double> phihat_var(nbetas, 0);
  std::vector<double> ncounts(nbetas, 0);
  std::vector<int> bin1(nbetas, 0);
  std::vector<int> bin2(nbetas, 0);
  
  cts_to_signal_mat_core(N, &cts_bin1[0], &cts_bin2[0], &count[0], &lmu_nosig[0], &weight[0], nbins, dispersion,
                         &phi[0], eCprime, &phihat[0], &phihat_var[0], &ncounts[0], &bin1[0], &bin2[0], diag_rm);
  
  IntegerVector bin1_i, bin2_i;
  NumericVector phihat_i, phihat_var_i, ncounts_i, weight_i, didx_i;
  bin1_i = wrap(bin1);
  bin2_i = wrap(bin2);
  bin1_i.attr("levels") = as<IntegerVector>(cts["bin1"]).attr("levels");
  bin2_i.attr("levels") = as<IntegerVector>(cts["bin2"]).attr("levels");
  bin1_i.attr("class") = CharacterVector::create("ordered", "factor");
  bin2_i.attr("class") = CharacterVector::create("ordered", "factor");
  phihat_i = wrap(phihat);
  phihat_var_i = wrap(phihat_var);
  ncounts_i = wrap(ncounts);
  weight_i = 1/phihat_var_i;
  didx_i = bin2_i-bin1_i;
  
  
  return DataFrame::create(_["bin1"]=bin1_i, _["bin2"]=bin2_i, _["phihat"]=phihat_i,
                           _["phihat.var"]=phihat_var_i, _["ncounts"]=ncounts_i, _["weight"]=weight_i,
                           _["diag.idx"]=didx_i);
}

std::vector<double> soft_threshold(const std::vector<double>& beta, double eCprime, double lam1) {
  std::vector<double> phi;
  phi.reserve(beta.size());
  for (std::vector<double>::const_iterator it = beta.begin(); it != beta.end(); ++it) {
    double val = *it - eCprime;
    phi.push_back((val > 0) ? std::max(0., val - lam1) : std::min(0., val + lam1));
  }
  return phi;
}

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                           double lam1, double lam2, double eCprime,
                           double alpha, double inflate, int ninner, double converge,
                           int diag_rm, NumericVector z_i, NumericVector u_i, NumericVector beta_i)
{
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<int> trails_r = as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
  std::vector<double> z_r = as<std::vector<double> >(z_i);
  std::vector<double> u_r = as<std::vector<double> >(u_i);
  std::vector<double> beta_r = as<std::vector<double> >(beta_i); //2d fused lasso before soft-thresholding
  std::vector<double> phi_r = soft_threshold(beta_r, eCprime, lam1); //sparse fused lasso soft-thresholds values
  std::vector<double> phi_old = phi_r;
  
  int step;
  int res=0;
  //printf(" Perf iteration: start with lam=%f alpha=%f phi[0]=%f z[0]=%f u[0]=%f\n",
  //       lam, alpha, phi_r[0], z_r[0], u_r[0]);
  for (step=1; step<=nouter; ++step) {
    //compute weights
    const DataFrame mat = cts_to_signal_mat(cts, nbins, dispersion, phi_r, eCprime, diag_rm);
    std::vector<double> y_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
    std::vector<double> w_r = Rcpp::as<std::vector<double> >(mat["weight"]);
    
    //compute fused lasso solution
    res += graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails, &trails_r[0], &breakpoints_r[0],
                                          lam2, &alpha, inflate, ninner, converge,
                                          &beta_r[0], &z_r[0], &u_r[0]);
    
    //soft-threshold it at the selected parameters
    phi_r = soft_threshold(beta_r, eCprime, lam1);
    
    //check convergence
    double maxval = std::abs(phi_r[0]-phi_old[0]);
    for (int i=1; i<N; ++i) maxval = std::max(std::abs(phi_r[i]-phi_old[i]), maxval);
    //printf(" Iteration %d with alpha=%f reached maxval=%.5e after %d steps (phi[0]=%f z[0]=%f u[0]=%f)\n",
    //       step,alpha,maxval,res,phi_r[0],z_r[0], u_r[0]);
    if (maxval<converge) break;
    phi_old = phi_r;
  }
  //printf(" Perf iteration: end   with lam=%f alpha=%f phi[0]=%f z[0]=%f u[0]=%f nouter=%d ninner=%d\n",
  //       lam, alpha, phi_r[0], z_r[0], u_r[0], step, res);
  return List::create(_["beta"]=wrap(beta_r), _["alpha"]=wrap(alpha),
                      _["mat"]=cts_to_signal_mat(cts, nbins, dispersion, phi_r, eCprime, diag_rm),
                      _["z"]=wrap(z_r), _["u"]=wrap(u_r), _["nouter"]=step, _["ninner"]=res);
}

List wgfl_signal_perf(const DataFrame cts, double dispersion, int nouter, int nbins,
                      int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                      double lam1, double lam2, double eCprime, double alpha, double inflate, int ninner,
                      double converge, int diag_rm)
{
  NumericVector z_i(breakpoints_i(ntrails-1));
  NumericVector u_i(breakpoints_i(ntrails-1));
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  NumericVector beta_i(N);
  //printf("Fused lasso cold perf iteration with %d coefficients\n",phi.size());
  return wgfl_signal_perf_warm(cts, dispersion, nouter, nbins, ntrails, trails_i, breakpoints_i,
                               lam1, lam2, eCprime, alpha, inflate, ninner, converge, diag_rm, z_i, u_i, beta_i);
}


