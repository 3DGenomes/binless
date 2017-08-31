#include <Rcpp.h>
#include <vector>

#include "DataLikelihoods.hpp"
#include "util.hpp"

Rcpp::NumericVector SignalLikelihood::get_chi_square(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> value_r = as<std::vector<double> >(value_);
        Rcpp::NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
        
        const Rcpp::NumericVector chisq = weight_ * SQUARE(valuehat_ - (soft + eCprime));
        return chisq;
    }

Rcpp::NumericVector DifferenceLikelihood::get_chi_square(double LB, double UB) const {
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    std::vector<double> value_r = as<std::vector<double> >(value_);
    std::vector<double> soft_r = soft_threshold(value_r, eCprime, lambda1);
    Rcpp::NumericVector soft = wrap(soft_r);
    std::vector<double> valuehat_r = as<std::vector<double> >(valuehat_);
    std::vector<double> valuehat_ref_r = as<std::vector<double> >(valuehat_ref_);
    std::vector<double> var, var_ref;
    var.reserve(valuehat_r.size());
    var_ref.reserve(valuehat_r.size());
    for (int i=0; i<valuehat_r.size(); ++i) {
        var.push_back(1/weight_(i));
        var_ref.push_back(1/weight_ref_(i));
    }
    Rcpp::NumericVector phi_ref = wrap(compute_phi_ref(soft_r, valuehat_r, var, valuehat_ref_r, var_ref));
    
    const Rcpp::NumericVector chisq = weight_ * SQUARE(valuehat_ - (soft + eCprime + phi_ref))
    + weight_ref_ * SQUARE(valuehat_ref_ - phi_ref);
    return chisq;
}
