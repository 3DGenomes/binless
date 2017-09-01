#ifndef LIKELIHOODS_HPP
#define LIKELIHOODS_HPP

#include <Rcpp.h>

#include "BinnedData.hpp"

class SignalLikelihood {
public:
    typedef Rcpp::NumericVector var_t;
    SignalLikelihood(const SignalBinnedData& data, const var_t& beta_phi) :
       beta_phi_(beta_phi), weight_(data.get_weight()), phihat_(data.get_phihat()) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const;
    
private:
    Rcpp::NumericVector beta_phi_, weight_, phihat_;
};

class DifferenceLikelihood {
public:
    typedef Rcpp::NumericVector var_t;
    DifferenceLikelihood(const DifferenceBinnedData& data, const var_t& beta_delta) :
    beta_delta_(beta_delta), weight_(data.get_weight()), deltahat_(data.get_deltahat()),
    weight_ref_(data.get_weight_ref()), phihat_ref_(data.get_phihat_ref()) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const;
    
private:
    Rcpp::NumericVector beta_delta_, weight_, deltahat_, weight_ref_, phihat_ref_;
};

#endif
