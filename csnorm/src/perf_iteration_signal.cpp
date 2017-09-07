#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "perf_iteration_signal.hpp"

#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "SignalWeightsUpdater.hpp"
#include "IRLSEstimator.hpp"
#include "CVEstimator.hpp"
#include "RawData.hpp"
#include "BinnedData.hpp"
#include "Traits.hpp"
#include "SparsityEstimator.hpp"

#include "util.hpp" //SQUARE
#include "cts_to_mat.hpp" //cts_to_signal_mat
#include "graph_helpers.hpp" //get_patch_numbers



List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                           double lam2, double alpha, double converge,
                           const List outliers, NumericVector beta_i) {
    //Class that holds all the data. Other classes reference to it.
    Signal::raw_t raw(nbins, dispersion, cts, outliers);
    Signal::binned_t binned; //stored here, but will be populated by WeightsUpdater
    //setup computation of fused lasso solution
    FusedLassoGaussianEstimator<GFLLibrary> flo(nbins, converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    SignalWeightsUpdater wt(raw,binned); //size of the problem and input data
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
    
    DataFrame finalmat = DataFrame::create(_["bin1"]=binned.get_bin1(),
                                           _["bin2"]=binned.get_bin2(),
                                           _["phihat"]=binned.get_phihat(),
                                           _["ncounts"]=binned.get_ncounts(),
                                           _["weight"]=binned.get_weight(),
                                           _["diag.idx"]=binned.get_diag_idx(),
                                           _["diag.grp"]=binned.get_diag_grp(),
                                           _["beta"]=beta,
                                           _["phi"]=beta);
    
    return List::create(_["beta"]=wrap(beta), _["alpha"]=wrap(alpha),
                        _["phi"]=wrap(beta), _["mat"]=finalmat,
                        _["nouter"]=step, _["ninner"]=res,
                        _["eCprime"]=0, _["lambda1"]=0);
}

List wgfl_signal_BIC(const DataFrame cts, double dispersion, int nouter, int nbins,
                     double lam2,  double alpha, double tol_val,
                     List outliers, NumericVector beta_i, double lambda1_min, int refine_num,
                     bool constrained, bool fixed) {
    
    //Class that holds all the data. Other classes reference to it.
    Signal::raw_t raw(nbins, dispersion, cts, outliers);
    Signal::binned_t binned; //stored here, but will be populated by WeightsUpdater
    //setup computation of fused lasso solution
    bool converged = true;
    const double converge = tol_val/20.;
    FusedLassoGaussianEstimator<GFLLibrary> flo(raw.get_nbins(), converge); //size of the problem and convergence criterion
    flo.setUp(alpha);
    SignalWeightsUpdater wt(raw,binned); //size of the problem and input data
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
    //retrieve statistics and compute patches
    alpha = flo.get_alpha();
    beta = flo.get();
    IntegerVector patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(),
                                              binned.get_bin2(), binned.get_beta_phi());
    binned.set_patchno(patchno);
    
    //compute CV datasets at optimized weights
    auto cv = make_CVEstimator(flo, wt, 0);
    cv.compute(beta, lam2);
    
    //optimize lambda1 and eC
    NumericVector opt;
    {
        NumericVector beta_cv = Rcpp::wrap(cv.get_beta_cv());
        IntegerVector cv_grp = Rcpp::wrap(cv.get_cvgroup());
        if (fixed) { // is eCprime fixed to 0?
            if (!constrained) stop("expected constrained==T when fixed==T");
            SparsityEstimator<Signal, CVkSD<1>, ZeroOffset, PositiveSign, ForbidDegeneracy> est(nbins, tol_val, binned, lam2, beta_cv, cv_grp);
            opt = est.optimize();
        } else {
            SparsityEstimator<Signal, CV, EstimatedOffset, PositiveSign, ForbidDegeneracy> est(nbins, tol_val, binned, lam2, beta_cv, cv_grp);
            opt = est.optimize();
        }
    }
    double lam1 = opt["lambda1"];
    double eCprime = opt["eCprime"];
    
    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(binned.get_beta());
    std::vector<double> phi_r = soft_threshold(beta_r, eCprime, lam1);

    //identify patches
    NumericVector phi = wrap(phi_r);
    patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(),
                                binned.get_bin2(), phi);

    //count the positive ones and deduce dof
    IntegerVector selected = patchno[abs(phi)>tol_val/2];
    const int dof = unique(selected).size();

    //retrieve BIC from previous computation
    const double BIC = opt["BIC"];
    const double BIC_sd = opt["BIC.sd"];
    
    DataFrame mat = DataFrame::create(_["bin1"]=binned.get_bin1(),
                                      _["bin2"]=binned.get_bin2(),
                                      _["phihat"]=binned.get_phihat(),
                                      _["ncounts"]=binned.get_ncounts(),
                                      _["weight"]=binned.get_weight(),
                                      _["diag.idx"]=binned.get_diag_idx(),
                                      _["diag.grp"]=binned.get_diag_grp(),
                                      _["beta"]=beta,
                                      _["phi"]=phi,
                                      _["beta_cv"]=cv.get_beta_cv(),
                                      _["cv.group"]=cv.get_cvgroup(),
                                      _["patchno"]=patchno);
    
    return List::create(_["phi"]=phi,
                        _["beta"]=beta_r, _["alpha"]=alpha, _["lambda2"]=lam2,
                        _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd, _["mat"]=mat, _["eCprime"]=eCprime,
                        _["lambda1"]=lam1, _["converged"]=converged);
}

