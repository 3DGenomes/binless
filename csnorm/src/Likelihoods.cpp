#include <Rcpp.h>
#include <vector>

#include "Likelihoods.hpp"
#include "util.hpp"

Rcpp::NumericVector SignalLikelihood::get_chi_square(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> beta_phi_r = Rcpp::as<std::vector<double> >(beta_phi_);
        Rcpp::NumericVector phi = Rcpp::wrap(soft_threshold(beta_phi_r, eCprime, lambda1));
        
        const Rcpp::NumericVector chisq = binned_.get_weight() * SQUARE(binned_.get_phihat() - (phi + eCprime));
        return chisq;
    }

Rcpp::NumericVector DifferenceLikelihood::get_chi_square(double LB, double UB) const {
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    std::vector<double> beta_delta_r = Rcpp::as<std::vector<double> >(beta_delta_);
    std::vector<double> delta_r = soft_threshold(beta_delta_r, eCprime, lambda1);
    Rcpp::NumericVector delta = Rcpp::wrap(delta_r);
    Rcpp::NumericVector phi_ref = compute_phi_ref(binned_, delta);
    const Rcpp::NumericVector chisq = binned_.get_weight() * SQUARE(binned_.get_phihat() - (delta + eCprime + phi_ref))
                                      + binned_.get_weight_ref() * SQUARE(binned_.get_phihat_ref() - phi_ref);
    return chisq;
}
