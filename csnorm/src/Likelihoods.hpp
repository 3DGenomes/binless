#ifndef LIKELIHOODS_HPP
#define LIKELIHOODS_HPP

#include <Rcpp.h>

#include "BinnedData.hpp"
#include "Traits.hpp"

//generic
template<typename Calculation> class Likelihood {};

//Signal
template<> class Likelihood<Signal> {
public:
    typedef Rcpp::NumericVector var_t;
    Likelihood(const BinnedData<Signal>& data, const var_t& beta_phi) : binned_(data), beta_phi_(beta_phi) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> beta_phi_r = Rcpp::as<std::vector<double> >(beta_phi_);
        Rcpp::NumericVector phi = Rcpp::wrap(soft_threshold(beta_phi_r, eCprime, lambda1));
        
        const Rcpp::NumericVector chisq = binned_.get_weight() * SQUARE(binned_.get_phihat() - (phi + eCprime));
        return chisq;
    }
    
private:
    const BinnedData<Signal>& binned_;
    const Rcpp::NumericVector beta_phi_;
};

//Difference
template<> class Likelihood<Difference> {
public:
    typedef Rcpp::NumericVector var_t;
    Likelihood(const BinnedData<Difference>& data, const var_t& beta_delta) :
    binned_(data), beta_delta_(beta_delta) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const {
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
    
private:
    const BinnedData<Difference>& binned_;
    const Rcpp::NumericVector beta_delta_;
};


#endif
