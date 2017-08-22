#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <ctime>

#include "perf_iteration_signal.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "GFLLibrary.hpp"

#include "util.hpp" //SQUARE
#include "cts_to_mat.hpp" //cts_to_signal_mat
#include "optimize_lambda1_eCprime.hpp" //cpp_optimize_lambda1_eCprime
#include "optimize_lambda1.hpp" //cpp_optimize_lambda1
#include "graph_helpers.hpp" //build_patch_graph_components


List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           double lam1, double lam2, double eCprime,
                           double alpha, double converge,
                           List outliers, NumericVector beta_i) {
    const int N = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<double> beta_r = as<std::vector<double> >
                                 (beta_i); //2d fused lasso before soft-thresholding
    std::vector<double> beta_old;
    std::vector<double> phi_r = soft_threshold(beta_r, eCprime,
                                lam1); //sparse fused lasso soft-thresholds values
    DataFrame mat;

    int step=0;
    double maxval=converge+1;
    std::clock_t c_start,c_end;
    double c_cts(0), c_gfl(0);
    /*Rcout << " Perf iteration: start with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
          << " lam1= " << lam1 << " eCprime= " << eCprime
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    /*Rcout << " eval init: lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/


    //setup computation of fused lasso solution, clamped at 50
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge);
    flo.setUp(alpha);
    
    while (step<=nouter & maxval>converge) {
        beta_old = beta_r;
        step++;
        //compute weights
        c_start = std::clock();
        mat = cts_to_signal_mat(cts, nbins, dispersion, phi_r, eCprime, outliers);
        std::vector<double> y_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
        std::vector<double> w_r = Rcpp::as<std::vector<double> >(mat["weight"]);
        c_end = std::clock();
        c_cts += c_end - c_start;

        //compute fused lasso
        c_start = std::clock();
        flo.optimize(y_r, beta_r, w_r, lam2);
        beta_r = flo.get();
        phi_r = flo.get(eCprime, lam1);
        alpha = flo.get_alpha();
        c_end = std::clock();
        c_gfl += c_end - c_start;

        //update residual
        maxval = std::abs(beta_r[0]-beta_old[0]);
        for (int i=1; i<N;
                ++i) maxval = std::max(std::abs(beta_r[i]-beta_old[i]), maxval);
        /*Rcout << " Iteration " << step << " with lam2= " << lam2 << " alpha= " << alpha << " reached maxval= " << maxval
              << " after " << flo.get_ninner() << " steps phi[0]= " << phi_r[0]
              << " lam1= " << lam1 << " eCprime= " << eCprime
              << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
        /*Rcout << " eval step "<< step << ": lam2= " << lam2 << " lam1= " << lam1 << " eCprime= " << eCprime
              << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
              << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r)))
              << " min(phihat)= " << min(NumericVector(wrap(y_r))) << " max(phihat)= "<< max(NumericVector(wrap(y_r))) << std::endl;*/


    }
    //if (step>nouter) Rcout << " warning: reached maximum number of outer iterations in wgfl_signal_perf_warm" << std::endl;
    /*Rcout << " Perf iteration: end   with lam2= " << lam2 << " alpha= " << alpha << " phi[0]= " << phi_r[0]
          << " lam1= " << lam1 << " eCprime= " << eCprime
          << " nouter= " << step << " ninner= " << res << " maxval= " << maxval
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    /*Rcout << " eval final: " << step << "/" << nouter << " lam2= " << lam2 << "lam1= " << lam1 << " eCprime= " << eCprime
          << " min(beta)= " << min(NumericVector(wrap(beta_r))) << " max(beta)= "<< max(NumericVector(wrap(beta_r)))
          << " min(phi)= " << min(NumericVector(wrap(phi_r))) << " max(phi)= "<< max(NumericVector(wrap(phi_r))) << std::endl;*/
    
    int res = flo.get_ninner();
    
    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["phihat"]=mat["phihat"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["weight"]=mat["weight"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta_r,
                                           _["phi"]=phi_r);
    
    return List::create(_["beta"]=wrap(beta_r), _["alpha"]=wrap(alpha),
                        _["phi"]=wrap(phi_r), _["mat"]=finalmat,
                        _["nouter"]=step, _["ninner"]=res,
                        _["eCprime"]=eCprime, _["lambda1"]=lam1, _["c_cts"]=c_cts, _["c_gfl"]=c_gfl);
}

List wgfl_signal_cv(const DataFrame mat, int nbins,
                    double lam2, double alpha, double converge, NumericVector beta_i) {
    const int N = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<double> beta_r = as<std::vector<double> >(beta_i);
    std::vector<double> phihat_r = as<std::vector<double> >(mat["phihat"]);
    std::vector<double> weight_r = as<std::vector<double> >(mat["weight"]);
    IntegerVector bin1 = as<IntegerVector>(mat["bin1"]);
    IntegerVector bin2 = as<IntegerVector>(mat["bin2"]);
    
    //build cv groups
    std::vector<int> cvgroup;
    const int ngroups=2;
    for (int i=0; i<N; ++i)
      cvgroup.push_back( (bin2[i]+bin1[i]) % ngroups ); // 2 cv groups in checkerboard pattern
    
    //setup computation of fused lasso solution, clamped at 50
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge);
    flo.setUp(alpha);
    
    //Compute fused lasso solutions on each group and report to beta_cv
    std::vector<double> beta_cv(N, -100);
    for (int g=0; g<ngroups; ++g) {
      //prepare data and weights for group g and copy initial values
      std::vector<double> p_r, w_r;
      for (int i=0; i<N; ++i) {
        if (cvgroup[i]==g) {
          p_r.push_back(0); //essential if lam2==0
          w_r.push_back(0);
        } else {
          p_r.push_back(phihat_r[i]);
          w_r.push_back(weight_r[i]);
        }
      }
      std::vector<double> values(beta_r);
      //compute fused lasso
      flo.optimize(p_r, values, w_r, lam2);
      values = flo.get();
      alpha = flo.get_alpha();
        
      //store fused solution at group positions back in beta_cv
      for (int i=0; i<N; ++i) if (cvgroup[i]==g) beta_cv[i] = values[i];
    }
    int res = flo.get_ninner();
    
    return List::create(_["beta_cv"]=wrap(beta_cv), _["cv.group"]=wrap(cvgroup),
                        _["ninner"]=wrap(res));
}

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int nbins,
                     double lam2,  double alpha, double tol_val,
                     List outliers, NumericVector beta_i, double lambda1_min, int refine_num,
                     bool constrained, bool fixed) {
    std::clock_t c_start,c_end;
    double c_cts(0), c_gfl(0), c_opt(0), c_init(0), c_brent(0), c_refine(0);
    double lam1=0, eCprime=0;
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_signal_perf_warm(cts, dispersion, nwarm, nbins, lam1, lam2, eCprime,
                                     alpha, tol_val/20., outliers, beta_i);
    c_cts += as<double>(ret["c_cts"]);
    c_gfl += as<double>(ret["c_gfl"]);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_signal_perf_warm(cts, dispersion, nouter, nbins, lam1, lam2, eCprime,
                                    alpha, tol_val/20., outliers, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
          //Rcout << " warning: cold start did not converge" <<std::endl;
          converged = false;
        }
        c_cts += as<double>(ret["c_cts"]);
        c_gfl += as<double>(ret["c_gfl"]);
    }
    
    //compute CV datasets at optimized weights
    DataFrame mat = as<DataFrame>(ret["mat"]);
    List cv_run = wgfl_signal_cv(mat, nbins, lam2,
                                 alpha, tol_val/20., beta_i);
      
    //optimize lambda1 and eC
    c_start = std::clock();
    DataFrame newmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["phihat"]=mat["phihat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["diag.idx"]=mat["diag.idx"],
                                         _["diag.grp"]=mat["diag.grp"],
                                         _["weight"]=mat["weight"],
                                         _["beta"]=ret["beta"],
                                         _["beta_cv"]=cv_run["beta_cv"],
                                         _["cv.group"]=cv_run["cv.group"],
                                         _["value"]=ret["beta"]);
    NumericVector opt;
    if (fixed) { // is eCprime fixed to 0?
        if (!constrained) stop("expected constrained==T when fixed==T");
        const bool positive = true;
        opt = cpp_optimize_lambda1(newmat, nbins, tol_val, positive, lambda1_min,
                                   refine_num);
    } else {
        opt = cpp_optimize_lambda1_eCprime(newmat, nbins, tol_val, constrained,
                                           lambda1_min, refine_num, lam2);
    }
    lam1 = opt["lambda1"];
    eCprime = opt["eCprime"];
    c_end = std::clock();

    c_opt += c_end - c_start;
    c_init += opt["c_init"];
    c_brent += opt["c_brent"];
    c_refine += opt["c_refine"];

    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(ret["beta"]);
    std::vector<double> phi_r = soft_threshold(beta_r, eCprime, lam1);

    //identify patches
    DataFrame submat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["valuehat"]=mat["phihat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["weight"]=mat["weight"],
                                         _["value"]=phi_r);
    List patches = build_patch_graph_components(nbins, submat, tol_val);

    //count the positive ones and deduce dof
    NumericVector phi = wrap(phi_r);
    IntegerVector patchno = patches["membership"];
    IntegerVector selected = patchno[abs(phi)>tol_val/2];
    const int dof = unique(selected).size();

    //retrieve BIC from previous computation
    const double BIC = opt["BIC"];
    const double BIC_sd = opt["BIC.sd"];
    
    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["phihat"]=mat["phihat"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["weight"]=mat["weight"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta_r,
                                           _["beta_cv"]=cv_run["beta_cv"],
                                           _["cv.group"]=cv_run["cv.group"],
                                           _["phi"]=phi_r,
                                           _["cv.group"]=cv_run["cv.group"],
                                           _["patchno"]=patchno);
    return List::create(_["phi"]=phi_r,
                        _["beta"]=beta_r, _["alpha"]=ret["alpha"], _["lambda2"]=lam2,
                        _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd, _["mat"]=finalmat, _["eCprime"]=eCprime,
                        _["lambda1"]=lam1,
                        _["c_cts"]=c_cts, _["c_gfl"]=c_gfl, _["c_opt"]=c_opt, _["c_init"]=c_init,
                        _["c_brent"]=c_brent, _["c_refine"]=c_refine, _["converged"]=converged);
}

List wgfl_signal_BIC_fixed(const DataFrame cts, double dispersion, int nouter, int nbins,
                     double lam1, double lam2, double eCprime, double alpha, double tol_val,
                     List outliers, NumericVector beta_i) {
    std::clock_t c_start,c_end;
    double c_cts(0), c_gfl(0), c_opt(0), c_init(0), c_brent(0), c_refine(0);
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_signal_perf_warm(cts, dispersion, nwarm, nbins, 0, lam2, 0,
                                     alpha, tol_val/20., outliers, beta_i);
    c_cts += as<double>(ret["c_cts"]);
    c_gfl += as<double>(ret["c_gfl"]);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_signal_perf_warm(cts, dispersion, nouter, nbins, lam1, lam2, eCprime,
                                    alpha, tol_val/20., outliers, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
            //Rcout << " warning: cold start did not converge" <<std::endl;
            converged=false;
        }
        c_cts += as<double>(ret["c_cts"]);
        c_gfl += as<double>(ret["c_gfl"]);
    }

    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(ret["beta"]);
    std::vector<double> phi_r = soft_threshold(beta_r, eCprime, lam1);

    //identify patches
    DataFrame mat = as<DataFrame>(ret["mat"]);
    DataFrame submat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["valuehat"]=mat["phihat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["weight"]=mat["weight"],
                                         _["value"]=phi_r);
    List patches = build_patch_graph_components(nbins, submat, tol_val);

    //count the positive ones and deduce dof
    NumericVector phi = wrap(phi_r);
    IntegerVector patchno = patches["membership"];
    IntegerVector selected = patchno[abs(phi)>tol_val/2];
    const int dof = unique(selected).size();

    //compute BIC
    NumericVector weight = submat["weight"];
    NumericVector phihat = submat["valuehat"];
    NumericVector ncounts = submat["ncounts"];
    const double BIC = sum(weight * SQUARE(phihat-(phi + eCprime))) + log(sum(
                           ncounts))*dof;

    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["phihat"]=mat["phihat"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["weight"]=mat["weight"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta_r,
                                           _["phi"]=phi_r,
                                           _["patchno"]=patchno);
    return List::create(_["phi"]=phi_r,
                        _["beta"]=beta_r, _["alpha"]=ret["alpha"], _["lambda2"]=lam2,
                        _["dof"]=dof, _["BIC"]=BIC, _["mat"]=finalmat, _["eCprime"]=eCprime,
                        _["lambda1"]=lam1,
                        _["c_cts"]=c_cts, _["c_gfl"]=c_gfl, _["converged"]=converged);
}

