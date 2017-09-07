#include <Rcpp.h>
#include <vector>

#include "Likelihoods.hpp"
#include "util.hpp"

Rcpp::NumericVector SignalLikelihood::get_chi_square(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> beta_phi_r = as<std::vector<double> >(beta_phi_);
        Rcpp::NumericVector phi = wrap(soft_threshold(beta_phi_r, eCprime, lambda1));
        
        const Rcpp::NumericVector chisq = binned_.get_weight() * SQUARE(binned_.get_phihat() - (phi + eCprime));
        return chisq;
    }

Rcpp::NumericVector DifferenceLikelihood::get_chi_square(double LB, double UB) const {
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    std::vector<double> beta_delta_r = as<std::vector<double> >(beta_delta_);
    std::vector<double> delta_r = soft_threshold(beta_delta_r, eCprime, lambda1);
    Rcpp::NumericVector delta = wrap(delta_r);
    std::vector<double> deltahat_r = as<std::vector<double> >(binned_.get_deltahat());
    Rcpp::NumericVector phihat_ref = binned_.get_phihat_ref();
    std::vector<double> phihat_ref_r = as<std::vector<double> >(phihat_ref);
    Rcpp::NumericVector weight = binned_.get_weight();
    Rcpp::NumericVector weight_ref = binned_.get_weight_ref();
    std::vector<double> phihat_r, var, var_ref;
    phihat_r.reserve(deltahat_r.size());
    var.reserve(deltahat_r.size());
    var_ref.reserve(deltahat_r.size());
    for (int i=0; i<deltahat_r.size(); ++i) {
        phihat_r.push_back(phihat_ref_r[i] + deltahat_r[i]);
        var.push_back(1/weight(i));
        var_ref.push_back(1/weight_ref(i));
    }
    Rcpp::NumericVector phi_ref = wrap(compute_phi_ref(delta_r, phihat_r, var, phihat_ref_r, var_ref));
    Rcpp::NumericVector phihat = wrap(phihat_r);
    const Rcpp::NumericVector chisq = weight * SQUARE(phihat - (delta + eCprime + phi_ref))
                                      + weight_ref * SQUARE(phihat_ref - phi_ref);
    return chisq;
}
