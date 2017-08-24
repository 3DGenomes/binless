#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "perf_iteration_diff.hpp"

#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "DifferenceWeightsUpdater.hpp"
#include "IRLSEstimator.hpp"

#include "cts_to_mat.hpp" //cts_to_diff_mat
#include "util.hpp" //SQUARE
#include "optimize_lambda1_diff.hpp" //cpp_optimize_lambda1
#include "graph_helpers.hpp" //build_patch_graph_components


List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins,
                         double lam2, double alpha, double converge,
                         List outliers, NumericVector phi_ref_i, NumericVector beta_i) {
    //setup computation of fused lasso solution
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    DifferenceWeightsUpdater wt(nbins, dispersion, cts, ref, outliers); //size of the problem and input data
    std::vector<double> beta = as<std::vector<double> >(beta_i);
    std::vector<double> phi_ref = as<std::vector<double> >(phi_ref_i);
    wt.setUp(phi_ref, beta); //initial guess of phi_ref provided here
                             //because we know the type of wt, but irls doesn't
    
    //do IRLS iterations until convergence
    auto irls = make_IRLSEstimator(converge, flo, wt); //number of iterations, convergence criterion and workers
    irls.optimize(nouter, beta, lam2);
    
    //retrieve statistics
    int res = flo.get_ninner();
    unsigned step = irls.get_nouter();
    alpha = flo.get_alpha();
    beta = flo.get();
    DataFrame mat = wt.get_mat();
    phi_ref = wt.get_phi_ref();
    
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
                                           _["beta"]=beta,
                                           _["delta"]=beta,
                                           _["phi_ref"]=phi_ref);
    
    return List::create(_["beta"]=wrap(beta), _["alpha"]=wrap(alpha),
                        _["phi.ref"]=wrap(phi_ref), _["delta"]=wrap(beta), _["mat"]=finalmat,
                        _["nouter"]=step, _["ninner"]=res,
                        _["eCprime"]=0, _["lambda1"]=0);
}

List wgfl_diff_cv(const DataFrame mat, int nbins,
                  double lam2, double alpha, double converge, NumericVector beta_i) {
    const int N = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<double> beta_r = as<std::vector<double> >(beta_i);
    std::vector<double> phihat_r = as<std::vector<double> >(mat["phihat"]);
    std::vector<double> phi_ref_r = as<std::vector<double> >(mat["phi_ref"]);
    std::vector<double> phihat_var_r = as<std::vector<double> >(mat["phihat.var"]);
    IntegerVector bin1 = as<IntegerVector>(mat["bin1"]);
    IntegerVector bin2 = as<IntegerVector>(mat["bin2"]);
    
    //build cv groups
    std::vector<int> cvgroup;
    const int ngroups=2;
    for (int i=0; i<N; ++i)
        cvgroup.push_back( (bin2[i]+bin1[i]) % ngroups ); // 2 cv groups in checkerboard pattern
    
    //setup computation of fused lasso solution, clamped at 50 (see Settings.hpp)
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge);
    flo.setUp(alpha);
    
    //Compute fused lasso solutions on each group and report to beta_cv
    std::vector<double> beta_cv(N, -100);
    for (int g=0; g<ngroups; ++g) {
        //prepare data and weights for group g and copy initial values
        std::vector<double> d_r, w_r;
        for (int i=0; i<N; ++i) {
            if (cvgroup[i]==g) {
                d_r.push_back(1000); //essential if lam2==0
                w_r.push_back(0);
            } else {
                d_r.push_back(phihat_r[i]-phi_ref_r[i]);
                w_r.push_back(1/phihat_var_r[i]);
            }
        }
        std::vector<double> values(beta_r);
        //compute fused lasso
        flo.optimize(d_r, values, w_r, lam2);
        values = flo.get();
        alpha = flo.get_alpha();
        
        //store fused solution at group positions back in beta_cv
        for (int i=0; i<N; ++i) if (cvgroup[i]==g) beta_cv[i] = values[i];
    }
    int res = flo.get_ninner();
    
    
    return List::create(_["beta_cv"]=wrap(beta_cv), _["cv.group"]=wrap(cvgroup),
                        _["ninner"]=wrap(res));
}



List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   double lam2,  double alpha, double tol_val,
                   List outliers, NumericVector phi_ref_i,  NumericVector beta_i, double lambda1_min,
                   int refine_num, bool constrained) {
    double lam1=0;
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_diff_perf_warm(cts, ref, dispersion, nwarm, nbins, lam2,
                                   alpha, tol_val/20., outliers, phi_ref_i, beta_i);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_diff_perf_warm(cts, ref, dispersion, nouter, nbins, lam2,
                                  alpha, tol_val/20., outliers, phi_ref_i, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
          //Rcout << " warning: cold start did not converge" <<std::endl;
          converged = false;
        }
    }
    
    //compute CV datasets at optimized weights
    DataFrame mat = as<DataFrame>(ret["mat"]);
    List cv_run = wgfl_diff_cv(mat, nbins, lam2, alpha, tol_val/20., beta_i);
    
    //optimize lambda1 assuming eCprime=0
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
    List patches = build_patch_graph_components(nbins, submat, tol_val);

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
    return List::create(_["phi.ref"]=phi_ref_r,
                        _["delta"]=delta, _["beta"]=beta_r,
                        _["alpha"]=ret["alpha"], _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd,
                        _["mat"]=finalmat, _["lambda1"]=lam1, _["eCprime"]=0, _["converged"]=converged);
}


List wgfl_diff_BIC_fixed(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   double lam1, double lam2, double alpha, double tol_val,
                   List outliers, NumericVector phi_ref_i,  NumericVector beta_i) {
    bool converged = true;
    //perf iteration for this set of values
    int nwarm = (int)(nouter/10.+1);
    List ret = wgfl_diff_perf_warm(cts, ref, dispersion, nwarm, nbins, lam2,
                                   alpha, tol_val/20., outliers, phi_ref_i, beta_i);
    //redo iteration if warm start did not work
    if (as<int>(ret["nouter"])>nwarm) {
        beta_i = NumericVector(beta_i.size(),0);
        //Rcout << " warning: warm start failed " << std::endl;
        ret = wgfl_diff_perf_warm(cts, ref, dispersion, nouter, nbins, lam2,
                                  alpha, tol_val/20., outliers, phi_ref_i, beta_i);
        if (as<int>(ret["nouter"])>nouter) {
          //Rcout << " warning: cold start did not converge" <<std::endl;
          converged = false;
        }
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
    List patches = build_patch_graph_components(nbins, submat, tol_val);

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
    return List::create(_["phi.ref"]=phi_ref,
                        _["delta"]=delta, _["beta"]=beta_r,
                        _["alpha"]=ret["alpha"], _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC,
                        _["mat"]=finalmat, _["lambda1"]=lam1, _["eCprime"]=0, _["converged"]=converged);
}





