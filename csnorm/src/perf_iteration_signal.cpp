#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "perf_iteration_signal.hpp"

#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "SignalWeightsUpdater.hpp"
#include "IRLSEstimator.hpp"

#include "util.hpp" //SQUARE
#include "cts_to_mat.hpp" //cts_to_signal_mat
#include "optimize_lambda1_eCprime.hpp" //cpp_optimize_lambda1_eCprime
#include "optimize_lambda1.hpp" //cpp_optimize_lambda1
#include "graph_helpers.hpp" //build_patch_graph_components


List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           double lam2, double alpha, double converge,
                           const List outliers, NumericVector beta_i) {
    //setup computation of fused lasso solution
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    SignalWeightsUpdater wt(nbins, dispersion, cts, outliers); //size of the problem and input data
    wt.setUp(); //for consistency. No-op, since there's no phi_ref to compute
    
    //do IRLS iterations until convergence
    auto irls = make_IRLSEstimator(converge, flo, wt); //number of iterations, convergence criterion and workers
    std::vector<double> beta = as<std::vector<double> >(beta_i);
    irls.optimize(nouter, beta, lam2);
    
    //retrieve statistics
    int res = flo.get_ninner();
    unsigned step = irls.get_nouter();
    alpha = flo.get_alpha();
    beta = flo.get();
    DataFrame mat = wt.get_mat();
    
    DataFrame finalmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["phihat"]=mat["phihat"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["weight"]=mat["weight"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta,
                                           _["phi"]=beta);
    
    return List::create(_["beta"]=wrap(beta), _["alpha"]=wrap(alpha),
                        _["phi"]=wrap(beta), _["mat"]=finalmat,
                        _["nouter"]=step, _["ninner"]=res,
                        _["eCprime"]=0, _["lambda1"]=0);
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
    
    //setup computation of fused lasso solution
    bool converged = true;
    const double converge = tol_val/20.;
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    SignalWeightsUpdater wt(nbins, dispersion, cts, outliers); //size of the problem and input data
    wt.setUp(); //for consistency. No-op, since there's no phi_ref to compute
    std::vector<double> beta = as<std::vector<double> >(beta_i);
    
    //do IRLS iterations until convergence
    //first, warm start
    unsigned nwarm = (unsigned)(nouter/10.+1);
    auto irls = make_IRLSEstimator(converge, flo, wt);
    irls.optimize(nwarm, beta, lam2);
    unsigned step = irls.get_nouter();
    //cold start if failed
    if (step>nwarm) {
        std::fill(beta.begin(), beta.end(), 0);
        irls.optimize(nouter, beta, lam2);
        step = irls.get_nouter();
        if (step>nouter) {
            //Rcout << " warning: cold start did not converge" <<std::endl;
            converged = false;
        }
    }
    //retrieve statistics
    alpha = flo.get_alpha();
    beta = flo.get();
    DataFrame mat = wt.get_mat();
    mat = DataFrame::create(_["bin1"]=mat["bin1"],
                                           _["bin2"]=mat["bin2"],
                                           _["phihat"]=mat["phihat"],
                                           _["phihat.var"]=mat["phihat.var"],
                                           _["ncounts"]=mat["ncounts"],
                                           _["weight"]=mat["weight"],
                                           _["diag.idx"]=mat["diag.idx"],
                                           _["diag.grp"]=mat["diag.grp"],
                                           _["beta"]=beta,
                                           _["phi"]=beta);
    
    //compute CV datasets at optimized weights
    List cv_run = wgfl_signal_cv(mat, nbins, lam2,
                                 alpha, tol_val/20., beta_i);
    
    //optimize lambda1 and eC
    DataFrame newmat = DataFrame::create(_["bin1"]=mat["bin1"],
                                         _["bin2"]=mat["bin2"],
                                         _["phihat"]=mat["phihat"],
                                         _["ncounts"]=mat["ncounts"],
                                         _["diag.idx"]=mat["diag.idx"],
                                         _["diag.grp"]=mat["diag.grp"],
                                         _["weight"]=mat["weight"],
                                         _["beta"]=mat["beta"],
                                         _["beta_cv"]=cv_run["beta_cv"],
                                         _["cv.group"]=cv_run["cv.group"],
                                         _["value"]=mat["beta"]);
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
    double lam1 = opt["lambda1"];
    double eCprime = opt["eCprime"];
    
    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(mat["beta"]);
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
                        _["beta"]=beta_r, _["alpha"]=alpha, _["lambda2"]=lam2,
                        _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd, _["mat"]=finalmat, _["eCprime"]=eCprime,
                        _["lambda1"]=lam1, _["converged"]=converged);
}

List wgfl_signal_BIC_fixed(const DataFrame cts, double dispersion, int nouter, int nbins,
                     double lam1, double lam2, double eCprime, double alpha, double tol_val,
                     List outliers, NumericVector beta_i) {
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_signal_perf_warm(cts, dispersion, nwarm, nbins, lam2,
                                     alpha, tol_val/20., outliers, beta_i);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_signal_perf_warm(cts, dispersion, nouter, nbins, lam2,
                                    alpha, tol_val/20., outliers, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
            //Rcout << " warning: cold start did not converge" <<std::endl;
            converged=false;
        }
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
                        _["lambda1"]=lam1, _["converged"]=converged);
}

