#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "perf_iteration_diff.hpp"

#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "DifferenceWeightsUpdater.hpp"
#include "IRLSEstimator.hpp"
#include "Preparation.hpp"
#include "RawData.hpp"
#include "BinnedData.hpp"
#include "Traits.hpp"
#include "SparsityEstimator.hpp"

#include "util.hpp" //SQUARE
#include "cts_to_mat.hpp" //cts_to_diff_mat
#include "graph_helpers.hpp" //get_patch_numbers


List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref,
                         double dispersion, int nouter, int nbins,
                         double lam2, double alpha, double converge,
                         List outliers, NumericVector phi_ref_i, NumericVector beta_i) {
    //Classes that hold all the data. Other classes reference to it.
    RawData<Difference> raw(nbins, dispersion, cts, ref, outliers);
    BinnedData<Difference> binned; //stored here, but will be populated by WeightsUpdater
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
    RawData<Difference> raw(nbins, dispersion, cts, ref, outliers);
    BinnedData<Difference> binned; //stored here, but will be populated by WeightsUpdater
    //setup computation of fused lasso solution
    bool converged = true;
    const double converge = tol_val/20.;
    FusedLassoGaussianEstimator<GFLLibrary> flo(raw.get_nbins(), converge); //size of the problem and convergence criterion
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
    //retrieve statistics and compute patches
    alpha = flo.get_alpha();
    beta = flo.get();
    IntegerVector patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(),
                                              binned.get_bin2(), binned.get_beta_delta());
    binned.set_patchno(patchno);
    
    //compute CV datasets at optimized weights
    auto cv = make_CVEstimator(flo, binned);
    cv.compute(beta, lam2);

    //optimize lambda1 assuming eCprime=0
    if (!constrained) stop("expected constrained==T when fixed==T");
    NumericVector opt;
    {
        NumericVector beta_cv = Rcpp::wrap(cv.get_beta_cv());
        IntegerVector cv_grp = Rcpp::wrap(cv.get_cvgroup());
        SparsityEstimator<Difference, CVkSD<1>, ZeroOffset, AnySign, ForbidDegeneracy> est(nbins, tol_val, binned, lam2, beta_cv, cv_grp);
        opt = est.optimize();
        
    }
    double lam1 = opt["lambda1"];
    
    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(binned.get_beta());
    std::vector<double> delta_r = soft_threshold(beta_r, 0, lam1);
    NumericVector delta = wrap(delta_r);
    
    //compute phi_ref
    Rcpp::NumericVector phi_ref_r = compute_phi_ref(binned, delta);
    
    //identify patches
    patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(),
                                binned.get_bin2(), delta);
    
    //count the positive ones and deduce dof
    IntegerVector selected = patchno[abs(delta)>tol_val/2];
    const int dof = unique(selected).size();

    //retrieve BIC from previous computation
    const double BIC = opt["BIC"];
    const double BIC_sd = opt["BIC.sd"];

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
                            _["delta"]=delta,
                            _["phi_ref"]=phi_ref_r,
                            _["beta_cv"]=cv.get_beta_cv(),
                            _["cv.group"]=cv.get_cvgroup(),
                            _["patchno"]=patchno);
    return List::create(_["phi.ref"]=phi_ref_r,
                        _["delta"]=delta, _["beta"]=beta,
                        _["alpha"]=alpha, _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd,
                        _["mat"]=mat, _["lambda1"]=lam1, _["eCprime"]=0, _["converged"]=converged);
}

