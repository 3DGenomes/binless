#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

#include "perf_iteration_diff.hpp"

#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "WeightsUpdater.hpp"
#include "IRLSEstimator.hpp"
#include "RawData.hpp"
#include "BinnedData.hpp"
#include "Traits.hpp"
#include "FixedSparsity.hpp"
#include "EstimatedSparsity.hpp"

#include "graph_helpers.hpp" //get_patch_numbers

List wgfl_diff_BIC(const DataFrame cts, const DataFrame ref, double dispersion,
                   int nouter, int nbins, List GFLState,
                   double lam2, double tol_val,
                   List metadata, NumericVector phi_ref_i,  NumericVector beta_i, bool constrained,
                   bool lambda1_fixed, double lambda1_fix_value) {
    
    //Classes that hold all the data. Other classes reference to it.
    RawData<Difference> raw(nbins, dispersion, cts, ref, metadata);
    BinnedData<Difference> binned; //stored here, but will be populated by WeightsUpdater
    //setup computation of fused lasso solution
    bool converged = true;
    const double converge = tol_val/20.;
    FusedLassoGaussianEstimator<GFLLibrary> flo(raw.get_nbins(), converge); //size of the problem and convergence criterion
    flo.set_state(GFLState);
    WeightsUpdater<Difference> wt(raw, binned); //size of the problem and input data
    std::vector<double> beta = as<std::vector<double> >(beta_i);
    std::vector<double> phi_ref = as<std::vector<double> >(phi_ref_i);
    wt.setUp(phi_ref, beta); //initial guess of phi_ref provided here
                             //because we know the type of wt, but irls doesn't
    
    //do IRLS iterations until convergence
    //first, warm start
    auto irls = make_IRLSEstimator(converge, flo, wt);
    irls.optimize(nouter, beta, lam2);
    converged = irls.has_converged();
    //cold start if failed
    if (!converged) {
        std::fill(beta.begin(), beta.end(), 0);
        flo.reset();
        irls.optimize(nouter, beta, lam2);
        converged = irls.has_converged();
        /*if (!converged) {
            Rcout << " warning: cold start did not converge" <<std::endl;
        }*/
    }
    
    //store GFL state for next round
    List new_GFLState = flo.get_state();
    
    //compute patches
    IntegerVector patchno = get_patch_numbers(nbins, tol_val, binned.get_bin1(),
                                              binned.get_bin2(), binned.get_beta_delta());
    binned.set_patchno(patchno);
    
    //optimize lambda1 assuming eCprime=0
    NumericVector opt;
    if (lambda1_fixed) {
        if (constrained) {
            auto est = make_FixedSparsity<BIC, ZeroOffset, AnySign, ForbidDegeneracy>(nbins, tol_val, binned, lam2, flo, lambda1_fix_value);
            opt = est.optimize();
        } else {
            auto est = make_FixedSparsity<BIC, ZeroOffset, AnySign, AllowDegeneracy>(nbins, tol_val, binned, lam2, flo, lambda1_fix_value);
            opt = est.optimize();
        }
    } else {
        if (constrained) {
          auto est = make_EstimatedSparsity<BIC, ZeroOffset, AnySign, ForbidDegeneracy>(nbins, tol_val, binned, lam2, flo);
          opt = est.optimize();
        } else {
          auto est = make_EstimatedSparsity<BIC, ZeroOffset, AnySign, AllowDegeneracy>(nbins, tol_val, binned, lam2, flo);
          opt = est.optimize();
        }
    }
    double lam1 = opt["lambda1"];
    
    //soft-threshold it at the selected parameters
    std::vector<double> beta_r = as<std::vector<double> >(binned.get_beta());
    std::vector<double> delta_r = soft_threshold(beta_r, 0, lam1);
    NumericVector delta = wrap(delta_r);
    
    //compute phi_ref
    Rcpp::NumericVector phi_ref_r = compute_phi_ref(binned, delta);
    
    //retrieve BIC from previous computation
    const double BIC = opt["BIC"];
    const double BIC_sd = opt["BIC.sd"];
    const unsigned dof = opt["dof"];

    DataFrame mat = DataFrame::create(_["bin1"]=binned.get_bin1(),
                            _["bin2"]=binned.get_bin2(),
                            _["phihat"]=binned.get_phihat(),
                            _["phihat.var"]=1/binned.get_weight(),
                            _["phihat.ref"]=binned.get_phihat_ref(),
                            _["phihat.var.ref"]=1/binned.get_weight_ref(),
                            _["nobs"]=binned.get_nobs(),
                            _["deltahat"]=binned.get_deltahat(),
                            _["weight"]=binned.get_weight(),
                            _["diag.idx"]=binned.get_diag_idx(),
                            _["diag.grp"]=binned.get_diag_grp(),
                            _["beta"]=binned.get_beta(),
                            _["delta"]=delta,
                            _["phi_ref"]=phi_ref_r,
                            _["patchno"]=binned.get_patchno());
    return List::create(_["phi.ref"]=phi_ref_r,
                        _["delta"]=delta, _["beta"]=binned.get_beta(), _["GFLState"]=new_GFLState,
                        _["lambda2"]=lam2, _["dof"]=dof, _["BIC"]=BIC, _["BIC.sd"]=BIC_sd,
                        _["mat"]=mat, _["lambda1"]=lam1, _["eCprime"]=0, _["converged"]=converged);
}

