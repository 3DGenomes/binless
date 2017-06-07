#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

#include "perf_iteration_signal.hpp"
#include "util.hpp" //SQUARE
#include "gfl_c.h" //cts_to_signal_mat_core
#include "graph_fl.h" //graph_fused_lasso_weight_warm
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include "optimize_lambda1_eCprime.hpp" //cpp_optimize_lambda1_eCprime


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

List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                           double lam1, double lam2, double eCprime,
                           double alpha, double inflate, int ninner, double converge,
                           int diag_rm, NumericVector beta_i)
{
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<int> trails_r = as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
  std::vector<double> beta_r = as<std::vector<double> >(beta_i); //2d fused lasso before soft-thresholding
  std::vector<double> beta_old;
  std::vector<double> phi_r = soft_threshold(beta_r, eCprime, lam1); //sparse fused lasso soft-thresholds values
  std::vector<double> u_r(trails_r.size(),0); //residuals set to zero
  std::vector<double> z_r;
  z_r.reserve(trails_r.size());
  for (int i=0; i<trails_r.size(); ++i) {z_r.push_back(beta_r[trails_r[i]]);} //z set to beta values along trails
  DataFrame mat;
  
  int step=0;
  int res=0;
  double maxval=converge+1;
  std::clock_t c_start,c_end;
  double c_cts(0), c_gfl(0);
  /*Rcout << " Perf iteration: start with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
        << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
  Rcout << " eval init: lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;
  
  
  while (step<=nouter & maxval>converge) {
    beta_old = beta_r;
    step++;
    //compute weights
    c_start = std::clock();
    mat = cts_to_signal_mat(cts, nbins, dispersion, phi_r, eCprime, diag_rm);
    std::vector<double> y_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
    std::vector<double> w_r = Rcpp::as<std::vector<double> >(mat["weight"]);
    c_end = std::clock();
    c_cts += c_end - c_start;
    
    //compute fused lasso solution
    c_start = std::clock();
    res += graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails, &trails_r[0], &breakpoints_r[0],
                                          lam2, &alpha, inflate, ninner, converge,
                                          &beta_r[0], &z_r[0], &u_r[0]);
    c_end = std::clock();
    c_gfl += c_end - c_start;
    
    //soft-threshold it at the selected parameters
    phi_r = soft_threshold(beta_r, eCprime, lam1);
    
    //update residual
    maxval = std::abs(beta_r[0]-beta_old[0]);
    for (int i=1; i<N; ++i) maxval = std::max(std::abs(beta_r[i]-beta_old[i]), maxval);
    /*Rcout << " Iteration " << step << " with lam2= " << lam2 << " alpha= " << alpha << " reached maxval= " << maxval
          << " after " << res << " steps phi[0]= " << phi_r[0]
          << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    /*Rcout << "eval step "<< step << ": lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
          << " min(phihat)= " << min(NumericVector(wrap(y_r))) << " max(beta)= "<< max(NumericVector(wrap(y_r)))
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    
    
  }
  if (step>nouter) Rcout << " warning: reached maximum number of outer iterations in wgfl_signal_perf_warm" << std::endl;
  /*Rcout << " Perf iteration: end   with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
        << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
        << " nouter= " << step << " ninner= " << res
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
  Rcout << " eval final: " << step << "/" << nouter << " lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;
  
  return List::create(_["beta"]=wrap(beta_r), _["alpha"]=wrap(alpha), _["phi"]=wrap(phi_r), _["mat"]=mat,
                      _["z"]=wrap(z_r), _["u"]=wrap(u_r), _["nouter"]=step, _["ninner"]=res,
                      _["eCprime"]=eCprime, _["lambda1"]=lam1, _["c_cts"]=c_cts, _["c_gfl"]=c_gfl);
}




List wgfl_signal_perf_opt_lambda1_eCprime(const DataFrame cts, double dispersion, int nouter, int opt_every,
                           int nbins, int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                           double lam2, double alpha, double inflate, int ninner, double converge,
                           int diag_rm,  double lam1_init, double eCprime_init, NumericVector beta_i,
                           double lambda1_min, int refine_num)
{
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  const bool constrained = true; //for signal step, constraint is always active
  double eCprime = eCprime_init, lam1 = lam1_init;
  double eCprime_old, lam1_old;
  std::vector<int> trails_r = as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
  std::vector<double> beta_r = as<std::vector<double> >(beta_i); //2d fused lasso before soft-thresholding
  std::vector<double> beta_old;
  std::vector<double> phi_r = soft_threshold(beta_r, eCprime, lam1); //sparse fused lasso soft-thresholds values
  std::vector<double> u_r(trails_r.size(),0); //residuals set to zero
  std::vector<double> z_r;
  z_r.reserve(trails_r.size());
  for (int i=0; i<trails_r.size(); ++i) {z_r.push_back(beta_r[trails_r[i]]);} //z set to beta values along trails
  DataFrame mat;
  
  int step=0;
  int res=0;
  double maxval=converge+1;
  std::clock_t c_start,c_end;
  double c_cts(0), c_gfl(0), c_opt(0), c_init(0), c_brent(0), c_refine(0);
  /*Rcout << " Perf iteration: start with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
        << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
  Rcout << "opt init: lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;
  while (step<=nouter & maxval>converge) {
    beta_old = beta_r;
    lam1_old = lam1;
    eCprime_old = eCprime;
    //update weights and lasso solution until convergence (or max steps reached), at fixed lambda1 and eCprime
    NumericVector beta_tmp = wrap(beta_r);
    List ret = wgfl_signal_perf_warm(cts, dispersion, opt_every, nbins, ntrails, trails_i, breakpoints_i, lam1, lam2, eCprime,
                                     alpha, inflate, ninner, converge, diag_rm, beta_i);
    int substep = ret["nouter"];
    beta_r = as<std::vector<double> >(ret["beta"]);
    phi_r = as<std::vector<double> >(ret["phi"]);
    mat = as<DataFrame>(ret["mat"]);
    c_cts += as<double>(ret["c_cts"]);
    c_gfl += as<double>(ret["c_gfl"]);
      
    //optimize eCprime and lambda1
    c_start = std::clock();
    DataFrame newmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["phihat"]=mat["phihat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["diag.idx"]=mat["diag.idx"],
                                         _["weight"]=mat["weight"],
                                         _["beta"]=beta_r,
                                         _["value"]=beta_r);
    NumericVector opt = cpp_optimize_lambda1_eCprime(newmat, nbins, converge*20, constrained, lambda1_min, refine_num);
    lam1 = opt["lambda1"];
    eCprime = opt["eCprime"];
    c_end = std::clock();
    
    c_opt += c_end - c_start;
    c_init += opt["c_init"];
    c_brent += opt["c_brent"];
    c_refine += opt["c_refine"];
    
    //soft-threshold it at the selected parameters
    phi_r = soft_threshold(beta_r, eCprime, lam1);
    
    //update residual
    maxval = std::max(std::abs(eCprime - eCprime_old), std::abs(lam1 - lam1_old));
    for (int i=0; i<N; ++i) maxval = std::max(std::abs(beta_r[i]-beta_old[i]), maxval);
    /*Rcout << " Iteration " << step << " with lam2= " << lam2 << " alpha= " << alpha << " reached maxval= " << maxval
          << " after " << res << " steps phi[0]= " << phi_r[0]
          << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    /*Rcout << "opt step "<< step << ": lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
          << " min(phihat)= " << min(as<NumericVector>(newmat["phihat"])) << " max(phihat)= "<< max(as<NumericVector>(newmat["phihat"]))
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    
    //update counter
    step += substep;
  }
  if (step>nouter) Rcout << " warning: reached maximum number of outer iterations in wgfl_signal_perf_opt_lambda1_eCprime " << std::endl;
  /*Rcout << " Perf iteration: end   with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
        << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
        << " nouter= " << step << " ninner= " << res
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << "max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
  Rcout << "opt final: " << step << "/" << nouter << " lam2= " << lam2 << " lam1= " << lam1 << " eCprime= " << eCprime
        << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
        << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;
  
  return List::create(_["beta"]=wrap(beta_r), _["alpha"]=wrap(alpha), _["phi"]=wrap(phi_r), _["mat"]=mat,
                      _["z"]=wrap(z_r), _["u"]=wrap(u_r), _["nouter"]=step, _["ninner"]=res,
                      _["eCprime"]=eCprime, _["lambda1"]=lam1, _["c_cts"]=c_cts, _["c_gfl"]=c_gfl, _["c_opt"]=c_opt,
                      _["c_init"]=c_init, _["c_brent"]=c_brent, _["c_refine"]=c_refine);
}

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int opt_every, int nbins,
                     int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                     double lam2,  double alpha, double inflate, int ninner, double tol_val,
                     int diag_rm, double lam1_init, double eCprime_init, NumericVector beta_i,
                     double lambda1_min, int refine_num) {
  
  //perf iteration for this set of values
  List ret = wgfl_signal_perf_opt_lambda1_eCprime(cts, dispersion, (int)(nouter/10.+1), opt_every, nbins, ntrails,
                                   trails_i, breakpoints_i,
                                   lam2, alpha, inflate, ninner, tol_val/20., diag_rm, lam1_init, eCprime_init,
                                   beta_i, lambda1_min, refine_num);
  //redo iteration if warm start did not work
  if (as<int>(ret["nouter"])>nouter) {
    Rcout << " warning: performing cold start due to failed warm start" <<std::endl;
    beta_i = NumericVector(beta_i.size(),0);
    ret = wgfl_signal_perf_opt_lambda1_eCprime(cts, dispersion, nouter, opt_every, nbins, ntrails, trails_i, breakpoints_i,
                                               lam2, alpha, inflate, ninner, tol_val/20., diag_rm, lambda1_min, 0,
                                               beta_i, lambda1_min, refine_num);
  }
  
  //identify patches
  DataFrame retmat = wrap(ret["mat"]);
  DataFrame submat = DataFrame::create(_["bin1"]=retmat["bin1"],
                                       _["bin2"]=retmat["bin2"],
                                       _["valuehat"]=retmat["phihat"],
                                       _["ncounts"]=retmat["ncounts"],
                                       _["weight"]=retmat["weight"],
                                       _["value"]=ret["phi"]);
  List patches = boost_build_patch_graph_components(nbins, submat, tol_val);
  
  //count the positive ones and deduce dof
  NumericVector phi = ret["phi"];
  IntegerVector patchno = patches["membership"];
  IntegerVector selected = patchno[abs(phi)>tol_val/2];
  const int dof = unique(selected).size();
  
  //compute BIC
  NumericVector weight = submat["weight"];
  NumericVector phihat = submat["valuehat"];
  NumericVector ncounts = submat["ncounts"];
  const double eCprime = ret["eCprime"];
  const double BIC = sum(weight * SQUARE(phihat-(phi + eCprime))) + log(sum(ncounts))*dof;
  
  DataFrame finalmat = DataFrame::create(_["bin1"]=retmat["bin1"],
                                       _["bin2"]=retmat["bin2"],
                                       _["valuehat"]=retmat["phihat"],
                                       _["ncounts"]=retmat["ncounts"],
                                       _["weight"]=retmat["weight"],
                                       _["value"]=ret["phi"],
                                       _["patchno"]=patchno);
  return List::create(_["z"]=ret["z"], _["u"]=ret["u"], _["beta"]=ret["beta"], _["alpha"]=ret["alpha"], _["lambda2"]=lam2,
                      _["dof"]=dof, _["BIC"]=BIC, _["mat"]=finalmat, _["eCprime"]=ret["eCprime"], _["lambda1"]=ret["lambda1"],
                      _["c_cts"]=ret["c_cts"], _["c_gfl"]=ret["c_gfl"], _["c_opt"]=ret["c_opt"],
                      _["c_init"]=ret["c_init"], _["c_brent"]=ret["c_brent"], _["c_refine"]=ret["c_refine"]);
}

