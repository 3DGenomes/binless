#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "perf_iteration_diff.hpp"

#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "DifferenceWeightsUpdater.hpp"
#include "IRLSEstimator.hpp"
#include "CVEstimator.hpp"
#include "RawData.hpp"
#include "BinnedData.hpp"
#include "Traits.hpp"
#include "Degeneracy.hpp"
#include "SparsityEstimator.hpp"

#include "cts_to_mat.hpp" //cts_to_diff_mat
#include "util.hpp" //SQUARE
#include "graph_helpers.hpp" //get_patch_numbers


List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins,
                         double lam2, double alpha, double converge,
                         List outliers, NumericVector phi_ref_i, NumericVector beta_i) {
    //Classes that hold all the data. Other classes reference to it.
    Difference::raw_t raw(nbins, dispersion, cts, ref, outliers);
    Difference::binned_t binned; //stored here, but will be populated by WeightsUpdater
    //setup computation of fused lasso solution
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    DifferenceWeightsUpdater wt(raw,binned); //size of the problem and input data
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
    phi_ref = wt.get_phi_ref();
    
    DataFrame mat = DataFrame::create(_["bin1"]=binned.get_bin1(),
                                      _["bin2"]=binned.get_bin2(),
                                      _["phihat"]=binned.get_phihat(),
                                      _["phihat.var"]=1/binned.get_weight(),
                                      _["phihat.ref"]=binned.get_phihat_ref(),
                                      _["phihat.var.ref"]=1/binned.get_weight_ref(),
                                      _["ncounts"]=binned.get_ncounts(),
                                      _["deltahat"]=binned.get_deltahat(),
                                      _["weight"]=binned.get_weight(),
                                      _["diag.idx"]=binned.get_diag_idx(),
                                      _["diag.grp"]=binned.get_diag_grp(),
                                      _["beta"]=beta,
                                      _["delta"]=beta,
                                      _["phi_ref"]=phi_ref);
    
    return List::create(_["beta"]=wrap(beta), _["alpha"]=wrap(alpha),
                        _["phi.ref"]=wrap(phi_ref), _["delta"]=wrap(beta), _["mat"]=mat,
                        _["nouter"]=step, _["ninner"]=res,
                        _["eCprime"]=0, _["lambda1"]=0);
}

List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins,
                   double lam2,  double alpha, double tol_val,
                   List outliers, NumericVector phi_ref_i,  NumericVector beta_i, double lambda1_min,
                   int refine_num, bool constrained) {
    
    //Classes that hold all the data. Other classes reference to it.
    Difference::raw_t raw(nbins, dispersion, cts, ref, outliers);
    Difference::binned_t binned; //stored here, but will be populated by WeightsUpdater
    //setup computation of fused lasso solution
    bool converged = true;
    const double converge = tol_val/20.;
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    DifferenceWeightsUpdater wt(raw, binned); //size of the problem and input data
    std::vector<double> beta = as<std::vector<double> >(beta_i);
    std::vector<double> phi_ref = as<std::vector<double> >(phi_ref_i);
    wt.setUp(phi_ref, beta); //initial guess of phi_ref provided here
                             //because we know the type of wt, but irls doesn't
    
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
    phi_ref = wt.get_phi_ref();
    
    
    //compute CV datasets at optimized weights
    auto cv = make_CVEstimator(flo, wt, 1000);
    cv.compute(beta, lam2);
    DataFrame mat = DataFrame::create(_["bin1"]=binned.get_bin1(),
                            _["bin2"]=binned.get_bin2(),
                            _["phihat"]=binned.get_phihat(),
                            _["phihat.var"]=1/binned.get_weight(),
                            _["phihat.ref"]=binned.get_phihat_ref(),
                            _["phihat.var.ref"]=1/binned.get_weight_ref(),
                            _["ncounts"]=binned.get_ncounts(),
                            _["deltahat"]=binned.get_deltahat(),
                            _["weight"]=binned.get_weight(),
                            _["diag.idx"]=binned.get_diag_idx(),
                            _["diag.grp"]=binned.get_diag_grp(),
                            _["beta"]=beta,
                            _["value"]=beta,
                            _["delta"]=beta,
                            _["phi_ref"]=phi_ref,
                            _["beta_cv"]=cv.get_beta_cv(),
                            _["cv.group"]=cv.get_cvgroup());
    
    //optimize lambda1 assuming eCprime=0
    if (!constrained) stop("expected constrained==T when fixed==T");
    NumericVector opt;
    {
        NumericVector beta_r = mat["beta"];
        NumericVector weight = mat["weight"];
        NumericVector phihat = mat["phihat"];
        NumericVector weight_ref = 1/as<NumericVector>(mat["phihat.var.ref"]);
        NumericVector phihat_ref = mat["phihat.ref"];
        NumericVector ncounts = mat["ncounts"];
        IntegerVector patchno = get_patch_numbers(nbins, mat, tol_val);
        NumericVector beta_cv = mat["beta_cv"];
        IntegerVector cv_grp = mat["cv.group"];
        DifferenceBinnedData data(beta_r, weight, phihat, weight_ref, phihat_ref, ncounts, patchno);
        SparsityEstimator<Difference, CVkSD<1>, ZeroOffset, AnySign, ForbidDegeneracy> est(nbins, tol_val, data, lam2, mat, beta_cv, cv_grp);
        opt = est.optimize();
        
    }
    double lam1 = opt["lambda1"];
    
    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(mat["beta"]);
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
    IntegerVector patchno = get_patch_numbers(nbins, submat, tol_val);

    //count the positive ones and deduce dof
    NumericVector delta = wrap(delta_r);
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
                                           _["beta_cv"]=mat["beta_cv"],
                                           _["cv.group"]=mat["cv.group"],
                                           _["delta"]=delta,
                                           _["phi.ref"]=phi_ref_r,
                                           _["patchno"]=patchno);
    return List::create(_["phi.ref"]=phi_ref_r,
                        _["delta"]=delta, _["beta"]=beta_r,
                        _["alpha"]=alpha, _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd,
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
    IntegerVector patchno = get_patch_numbers(nbins, submat, tol_val);

    //count the positive ones and deduce dof
    NumericVector delta = wrap(delta_r);
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





