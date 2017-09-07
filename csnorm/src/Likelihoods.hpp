#ifndef LIKELIHOODS_HPP
#define LIKELIHOODS_HPP

#include <Rcpp.h>

#include "BinnedData.hpp"

class SignalLikelihood {
public:
    typedef Rcpp::NumericVector var_t;
    SignalLikelihood(const SignalBinnedData& data, const var_t& beta_phi) : binned_(data), beta_phi_(beta_phi) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const;
    
private:
    const SignalBinnedData& binned_;
    const Rcpp::NumericVector beta_phi_;
};

class DifferenceLikelihood {
public:
    typedef Rcpp::NumericVector var_t;
    DifferenceLikelihood(const DifferenceBinnedData& data, const var_t& beta_delta) :
    binned_(data), beta_delta_(beta_delta) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const;
    
private:
    const DifferenceBinnedData& binned_;
    const Rcpp::NumericVector beta_delta_;
};

#endif
