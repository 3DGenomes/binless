#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

#include "perf_iteration_diff.hpp"
#include "perf_iteration_signal.hpp" //cts_to_signal_mat
#include "util.hpp" //SQUARE
#include "graph_fl.h" //graph_fused_lasso_weight_warm
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include "optimize_lambda1_diff.hpp" //cpp_optimize_lambda1



DataFrame cts_to_diff_mat(const DataFrame cts, const DataFrame ref, int nbins,
                          double dispersion,
                          std::vector<double>& phi_ref, std::vector<double>& delta, int diag_rm) {
    //assume eCprime = 0 for difference step
    const double eCprime=0;
    const DataFrame mat_ref = cts_to_signal_mat(ref, nbins, dispersion, phi_ref,
                              eCprime, diag_rm);

    std::vector<double> phi_oth;
    phi_oth.reserve(delta.size());
    for (int i=0; i<delta.size(); ++i) phi_oth[i] = phi_ref[i] + delta[i];
    const DataFrame mat_oth = cts_to_signal_mat(cts, nbins, dispersion, phi_oth,
                              eCprime, diag_rm);
    IntegerVector bin1 = mat_ref["bin1"];
    IntegerVector bin2 = mat_ref["bin2"];
    NumericVector phihat_ref = mat_ref["phihat"];
    NumericVector phihat_var_ref = mat_ref["phihat.var"];
    NumericVector phihat = mat_oth["phihat"];
    NumericVector phihat_var = mat_oth["phihat.var"];
    NumericVector ncounts_oth = mat_oth["ncounts"];
    NumericVector ncounts_ref = mat_ref["ncounts"];
    NumericVector deltahat = phihat-phihat_ref;
    NumericVector deltahat_var = phihat_var+phihat_var_ref;
    NumericVector ncounts = ncounts_ref+ncounts_oth;
    NumericVector weight = 1/deltahat_var;
    IntegerVector didx = mat_ref["diag.idx"];
    IntegerVector dgrp = mat_ref["diag.grp"];
    return DataFrame::create(_["bin1"]=bin1, _["bin2"]=bin2,
                             _["phihat"]=phihat, _["phihat.var"]=phihat_var,
                             _["phihat.ref"]=phihat_ref, _["phihat.var.ref"]=phihat_var_ref,
                             _["deltahat"]=deltahat, _["deltahat.var"]=deltahat_var, _["ncounts"]=ncounts,
                             _["weight"]=weight, _["diag.idx"]=didx, _["diag.grp"]=dgrp);
}

List wgfl_diff_cv(const DataFrame mat, int nbins,
                    int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                    double lam2, double alpha, double inflate, int ninner, double converge, NumericVector beta_i) {
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<int> trails_r = as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
  std::vector<double> beta_r = as<std::vector<double> >(beta_i);
  std::vector<double> deltahat_r = as<std::vector<double> >(mat["deltahat"]);
  std::vector<double> weight_r = as<std::vector<double> >(mat["weight"]);
  IntegerVector bin1 = as<IntegerVector>(mat["bin1"]);
  IntegerVector bin2 = as<IntegerVector>(mat["bin2"]);
  
  std::vector<double> u_r(trails_r.size(),0); //residuals set to zero
  std::vector<double> z_r;
  z_r.reserve(trails_r.size());
  for (int i=0; i<trails_r.size(); ++i) {
    z_r.push_back(beta_r[trails_r[i]]);   //z set to beta values along trails
  }
  
  //build cv groups
  std::vector<int> cvgroup;
  const int ngroups=2;
  for (int i=0; i<N; ++i)
    cvgroup.push_back( (bin2[i]+bin1[i]) % ngroups ); // 2 cv groups in checkerboard pattern
  
  //Compute fused lasso solutions on each group and report to beta_cv
  int res=0;
  std::vector<double> beta_cv(N, -100);
  for (int g=0; g<ngroups; ++g) {
    //prepare data and weights for group g and copy initial values
    std::vector<double> d_r, w_r;
    for (int i=0; i<N; ++i) {
      if (cvgroup[i]==g) {
        d_r.push_back(0); //essential if lam2==0
        w_r.push_back(0);
      } else {
        d_r.push_back(deltahat_r[i]);
        w_r.push_back(weight_r[i]);
      }
    }
    std::vector<double> values(beta_r);
    //compute fused lasso solution
    res += graph_fused_lasso_weight_warm (N, &d_r[0], &w_r[0], ntrails,
                                          &trails_r[0], &breakpoints_r[0],
                                                                      lam2, &alpha, inflate, ninner, converge,
                                                                      &values[0], &z_r[0], &u_r[0]);
    //store fused solution at group positions back in beta_cv
    for (int i=0; i<N; ++i) if (cvgroup[i]==g) beta_cv[i] = values[i];
  }
  
  return List::create(_["beta_cv"]=wrap(beta_cv), _["cv.group"]=wrap(cvgroup),
                      _["ninner"]=wrap(res));
}


List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins,
                         int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                         double lam1, double lam2, double alpha, double inflate, int ninner,
                         double converge,
                         int diag_rm, NumericVector phi_ref_i, NumericVector beta_i) {
    const int N = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<int> trails_r = as<std::vector<int> >(trails_i);
    std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
    std::vector<double> phi_ref_r = as<std::vector<double> >(phi_ref_i);
    std::vector<double> beta_r = as<std::vector<double> >(beta_i);
    std::vector<double> beta_old = beta_r;
    std::vector<double> delta_r = soft_threshold(beta_r, 0, lam1);
    std::vector<double> u_r(trails_r.size(),0); //residuals set to zero
    std::vector<double> z_r;
    z_r.reserve(trails_r.size());
    for (int i=0; i<trails_r.size(); ++i) {
        z_r.push_back(beta_r[trails_r[i]]);   //z set to beta values along trails
    }
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
    /*Rcout << " eval init: lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/


    while (step<=nouter & maxval>converge) {
        beta_old = beta_r;
        step++;
        //compute weights
        c_start = std::clock();
        mat = cts_to_diff_mat(cts, ref, nbins, dispersion, phi_ref_r, delta_r, diag_rm);
        std::vector<double> phihat_ref = Rcpp::as<std::vector<double> >
                                         (mat["phihat.ref"]);
        std::vector<double> phihat_var_ref = Rcpp::as<std::vector<double> >
                                             (mat["phihat.var.ref"]);
        std::vector<double> phihat = Rcpp::as<std::vector<double> >(mat["phihat"]);
        std::vector<double> phihat_var = Rcpp::as<std::vector<double> >
                                         (mat["phihat.var"]);
        std::vector<double> y_r, w_r;
        y_r.reserve(N);
        w_r.reserve(N);
        for (int i=0; i<N; ++i) {
            y_r.push_back(phihat[i]-phi_ref_r[i]);
            w_r.push_back(1/phihat_var[i]);
        }
        c_end = std::clock();
        c_cts += c_end - c_start;

        //compute fused lasso solution
        c_start = std::clock();
        res += graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails,
                                              &trails_r[0], &breakpoints_r[0],
                                              lam2, &alpha, inflate, ninner, converge,
                                              &beta_r[0], &z_r[0], &u_r[0]);
        const double beta_max = 50;
        for (std::vector<double>::iterator it = beta_r.begin(); it != beta_r.end(); ++it)
          *it = std::min(beta_max, std::max(-beta_max, *it));
        c_end = std::clock();
        c_gfl += c_end - c_start;

        //soft-threshold it at the selected parameters
        delta_r = soft_threshold(beta_r, 0, lam1);

        //update phi_ref
        phi_ref_r = compute_phi_ref(delta_r, phihat, phihat_var, phihat_ref,
                                    phihat_var_ref);
        
        //update residual
        maxval = std::abs(beta_r[0]-beta_old[0]);
        for (int i=1; i<N; ++i)
            maxval = std::max(std::abs(beta_r[i]-beta_old[i]), maxval);
        /*Rcout << " Iteration " << step << " with lam2= " << lam2 << " alpha= " << alpha << " reached maxval= " << maxval
              << " after " << res << " steps phi[0]= " << phi_r[0]
              << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
              << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
        /*Rcout << " eval step "<< step << ": lam2= " << lam2 << " lam1= " << lam1 << " eCprime= " << 0
              << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
              << " min(phi_ref)= " << min(NumericVector(wrap(phi_ref_r))) << " max(phi_ref)= "<< max(NumericVector(wrap(phi_ref_r)))
              << " min(phi)= " << min(NumericVector(wrap(delta_r))) << " max(delta)= "<< max(NumericVector(wrap(delta_r)))
              << " min(phihat)= " << min(NumericVector(wrap(y_r))) << " max(phihat)= "<< max(NumericVector(wrap(y_r))) << std::endl;*/


    }
    //if (step>nouter) Rcout << " warning: reached maximum number of outer iterations in wgfl_signal_perf_warm" << std::endl;
    /*Rcout << " Perf iteration: end   with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
          << " z[0]= " << z_r[0] << " u[0]= " << u_r[0] << " lam1= " << lam1 << " eCprime= " << eCprime
          << " nouter= " << step << " ninner= " << res
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    /*Rcout << " eval final: " << step << "/" << nouter << " lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/

    
    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["phihat"]=mat["phihat"],
                                           _["phihat.var"]=mat["phihat.var"],
                                           _["phihat.ref"]=mat["phihat.ref"],
                                           _["phihat.var.ref"]=mat["phihat.var.ref"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["deltahat"]=mat["deltahat"],
                                           _["weight"]=mat["weight"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta_r,
                                           _["delta"]=delta_r,
                                           _["phi_ref"]=phi_ref_r);
    
    return List::create(_["beta"]=wrap(beta_r), _["alpha"]=wrap(alpha),
                        _["phi.ref"]=wrap(phi_ref_r), _["delta"]=wrap(delta_r), _["mat"]=finalmat,
                        _["z"]=wrap(z_r), _["u"]=wrap(u_r), _["nouter"]=step, _["ninner"]=res,
                        _["eCprime"]=0, _["lambda1"]=lam1, _["c_cts"]=c_cts, _["c_gfl"]=c_gfl);
}


List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                   double lam2,  double alpha, double inflate, int ninner, double tol_val,
                   int diag_rm, NumericVector phi_ref_i,  NumericVector beta_i, double lambda1_min,
                   int refine_num, bool constrained) {
    std::clock_t c_start,c_end;
    double c_cts(0), c_gfl(0), c_opt(0), c_init(0), c_brent(0), c_refine(0);
    double lam1=0;
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_diff_perf_warm(cts, ref, dispersion, nwarm, nbins, ntrails,
                                   trails_i, breakpoints_i, lam1, lam2,
                                   alpha, inflate, ninner, tol_val/20., diag_rm, phi_ref_i, beta_i);
    c_cts += as<double>(ret["c_cts"]);
    c_gfl += as<double>(ret["c_gfl"]);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_diff_perf_warm(cts, ref, dispersion, nouter, nbins, ntrails,
                                  trails_i, breakpoints_i, lam1, lam2,
                                  alpha, inflate, ninner, tol_val/20., diag_rm, phi_ref_i, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
          //Rcout << " warning: cold start did not converge" <<std::endl;
          converged = false;
        }
        c_cts += as<double>(ret["c_cts"]);
        c_gfl += as<double>(ret["c_gfl"]);
    }
    
    //compute CV datasets at optimized weights
    DataFrame mat = as<DataFrame>(ret["mat"]);
    List cv_run = wgfl_diff_cv(mat, nbins, ntrails, trails_i, breakpoints_i, lam2,
                               alpha, inflate, ninner, tol_val/20., beta_i);
    
    //optimize lambda1 assuming eCprime=0
    c_start = std::clock();
    DataFrame newmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["phihat"]=mat["phihat"],
                                         _["phihat.ref"]=mat["phihat.ref"],
                                         _["phihat.var"]=mat["phihat.var"],
                                         _["phihat.var.ref"]=mat["phihat.var.ref"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["diag.idx"]=mat["diag.idx"],
                                         _["diag.grp"]=mat["diag.grp"],
                                         _["weight"]=mat["weight"],
                                         _["beta"]=ret["beta"],
                                         _["value"]=ret["beta"],
                                         _["beta_cv"]=cv_run["beta_cv"],
                                         _["cv.group"]=cv_run["cv.group"]);
    if (!constrained) stop("expected constrained==T when fixed==T");
    NumericVector opt = cpp_optimize_lambda1_diff(newmat, nbins, tol_val,
                        lambda1_min, refine_num);
    lam1 = opt["lambda1"];
    c_end = std::clock();
    c_opt += c_end - c_start;
    c_init += opt["c_init"];
    c_brent += opt["c_brent"];
    c_refine += opt["c_refine"];

    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(ret["beta"]);
    std::vector<double> delta_r = soft_threshold(beta_r, 0, lam1);

    //compute phi_ref
    std::vector<double> phihat_ref_r = Rcpp::as<std::vector<double> >
                                       (mat["phihat.ref"]);
    std::vector<double> phihat_var_ref_r = Rcpp::as<std::vector<double> >
                                           (mat["phihat.var.ref"]);
    std::vector<double> phihat_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
    std::vector<double> phihat_var_r = Rcpp::as<std::vector<double> >
                                       (mat["phihat.var"]);
    std::vector<double> phi_ref_r = compute_phi_ref(delta_r, phihat_r, phihat_var_r,
                                    phihat_ref_r, phihat_var_ref_r);

    //identify patches
    DataFrame submat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["valuehat"]=mat["deltahat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["weight"]=mat["weight"],
                                         _["value"]=delta_r);
    List patches = boost_build_patch_graph_components(nbins, submat, tol_val);

    //count the positive ones and deduce dof
    NumericVector delta = wrap(delta_r);
    IntegerVector patchno = patches["membership"];
    IntegerVector selected = patchno[abs(delta)>tol_val/2];
    const int dof = unique(selected).size();

    //retrieve BIC from previous computation
    const double BIC = opt["BIC"];
    const double BIC_sd = opt["BIC.sd"];

    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["deltahat"]=mat["deltahat"],
                                           _["weight"]=mat["weight"],
                                           _["phihat.ref"]=mat["phihat.ref"],
                                           _["phihat.var.ref"]=mat["phihat.var.ref"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta_r,
                                           _["beta_cv"]=cv_run["beta_cv"],
                                           _["cv.group"]=cv_run["cv.group"],
                                           _["delta"]=delta,
                                           _["phi.ref"]=phi_ref_r,
                                           _["patchno"]=patchno);
    return List::create(_["z"]=ret["z"], _["u"]=ret["u"], _["phi.ref"]=phi_ref_r,
                        _["delta"]=delta, _["beta"]=beta_r,
                        _["alpha"]=ret["alpha"], _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd,
                        _["mat"]=finalmat, _["lambda1"]=lam1, _["eCprime"]=0,
                        _["c_cts"]=c_cts, _["c_gfl"]=c_gfl, _["c_opt"]=c_opt, _["c_init"]=c_init,
                        _["c_brent"]=c_brent, _["c_refine"]=c_refine, _["converged"]=converged);
}


List wgfl_diff_BIC_fixed(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                   double lam1, double lam2, double alpha, double inflate, int ninner, double tol_val,
                   int diag_rm, NumericVector phi_ref_i,  NumericVector beta_i) {
    std::clock_t c_start,c_end;
    double c_cts(0), c_gfl(0), c_opt(0), c_init(0), c_brent(0), c_refine(0);
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_diff_perf_warm(cts, ref, dispersion, nwarm, nbins, ntrails,
                                   trails_i, breakpoints_i, 0, lam2,
                                   alpha, inflate, ninner, tol_val/20., diag_rm, phi_ref_i, beta_i);
    c_cts += as<double>(ret["c_cts"]);
    c_gfl += as<double>(ret["c_gfl"]);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_diff_perf_warm(cts, ref, dispersion, nouter, nbins, ntrails,
                                  trails_i, breakpoints_i, 0, lam2,
                                  alpha, inflate, ninner, tol_val/20., diag_rm, phi_ref_i, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
          //Rcout << " warning: cold start did not converge" <<std::endl;
          converged = false;
        }
        c_cts += as<double>(ret["c_cts"]);
        c_gfl += as<double>(ret["c_gfl"]);
    }

    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(ret["beta"]);
    std::vector<double> delta_r = soft_threshold(beta_r, 0, lam1);

    //compute phi_ref
    DataFrame mat = as<DataFrame>(ret["mat"]);
    std::vector<double> phihat_ref_r = Rcpp::as<std::vector<double> >
                                       (mat["phihat.ref"]);
    std::vector<double> phihat_var_ref_r = Rcpp::as<std::vector<double> >
                                           (mat["phihat.var.ref"]);
    std::vector<double> phihat_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
    std::vector<double> phihat_var_r = Rcpp::as<std::vector<double> >
                                       (mat["phihat.var"]);
    std::vector<double> phi_ref_r = compute_phi_ref(delta_r, phihat_r, phihat_var_r,
                                    phihat_ref_r, phihat_var_ref_r);

    //identify patches
    DataFrame submat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["valuehat"]=mat["deltahat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["weight"]=mat["weight"],
                                         _["value"]=delta_r);
    List patches = boost_build_patch_graph_components(nbins, submat, tol_val);

    //count the positive ones and deduce dof
    NumericVector delta = wrap(delta_r);
    IntegerVector patchno = patches["membership"];
    IntegerVector selected = patchno[abs(delta)>tol_val/2];
    const int dof = unique(selected).size();

    //compute BIC
    NumericVector weight = mat["weight"];
    NumericVector deltahat = mat["deltahat"];
    NumericVector phihat_ref = mat["phihat.ref"];
    NumericVector phi_ref = wrap(phi_ref_r);
    NumericVector phihat_var_ref = mat["phihat.var.ref"];
    NumericVector ncounts = mat["ncounts"];
    const double BIC = sum(weight * SQUARE(deltahat - delta) + SQUARE(
                               phihat_ref - phi_ref)/phihat_var_ref) + log(sum(ncounts))*dof;

    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["deltahat"]=mat["deltahat"],
                                           _["weight"]=mat["weight"],
                                           _["phihat.ref"]=mat["phihat.ref"],
                                           _["phihat.var.ref"]=mat["phihat.var.ref"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta_r,
                                           _["delta"]=delta,
                                           _["phi.ref"]=phi_ref,
                                           _["patchno"]=patchno);
    return List::create(_["z"]=ret["z"], _["u"]=ret["u"], _["phi.ref"]=phi_ref,
                        _["delta"]=delta, _["beta"]=beta_r,
                        _["alpha"]=ret["alpha"], _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC,
                        _["mat"]=finalmat, _["lambda1"]=lam1, _["eCprime"]=0,
                        _["c_cts"]=c_cts, _["c_gfl"]=c_gfl, _["converged"]=converged);
}





