#ifndef LIKELIHOODS_HPP
#define LIKELIHOODS_HPP

#include <Rcpp.h>

#include "BinnedData.hpp"

class SignalLikelihood {
public:
    typedef Rcpp::NumericVector var_t;
    SignalLikelihood(const SignalBinnedData& data, const var_t& value) :
       value_(value), weight_(data.get_weight()), valuehat_(data.get_valuehat()) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const;
    
private:
    Rcpp::NumericVector value_, weight_, valuehat_;
};

class DifferenceLikelihood {
public:
    typedef Rcpp::NumericVector var_t;
    DifferenceLikelihood(const DifferenceBinnedData& data, const var_t& value) :
    value_(value), weight_(data.get_weight()), valuehat_(data.get_valuehat()),
    weight_ref_(data.get_weight_ref()), valuehat_ref_(data.get_valuehat_ref()) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const;
    
private:
    Rcpp::NumericVector value_, weight_, valuehat_, weight_ref_, valuehat_ref_;
};

#endif
