#ifndef MODELS_HPP
#define MODELS_HPP

#include "util.hpp"
#include "Data.hpp"

class DifferenceLikelihood {
public:
    DifferenceLikelihood(const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
                         const NumericVector& weight_ref, const NumericVector& valuehat_ref) :
       value_(value), weight_(weight), valuehat_(valuehat), weight_ref_(weight_ref),
       valuehat_ref_(valuehat_ref) {}
    
    DifferenceLikelihood(const DifferenceData& data) :
       value_(data.get_value()), weight_(data.get_weight()), valuehat_(data.get_valuehat()),
       weight_ref_(data.get_weight_ref()), valuehat_ref_(data.get_valuehat_ref()) {}
    
    NumericVector get_chi_square(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> value_r = as<std::vector<double> >(value_);
        std::vector<double> soft_r = soft_threshold(value_r, eCprime, lambda1);
        NumericVector soft = wrap(soft_r);
        std::vector<double> valuehat_r = as<std::vector<double> >(valuehat_);
        std::vector<double> valuehat_ref_r = as<std::vector<double> >(valuehat_ref_);
        std::vector<double> var, var_ref;
        var.reserve(valuehat_r.size());
        var_ref.reserve(valuehat_r.size());
        for (int i=0; i<valuehat_r.size(); ++i) {
            var.push_back(1/weight_(i));
            var_ref.push_back(1/weight_ref_(i));
        }
        NumericVector phi_ref = wrap(compute_phi_ref(soft_r, valuehat_r, var, valuehat_ref_r, var_ref));
        
        const NumericVector chisq = weight_ * SQUARE(valuehat_ - (soft + eCprime + phi_ref))
                                       + weight_ref_ * SQUARE(valuehat_ref_ - phi_ref);
        return chisq;
    }

private:
    NumericVector value_, weight_, valuehat_, weight_ref_, valuehat_ref_;
};

class SignalLikelihood {
public:
    SignalLikelihood(const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat) :
       value_(value), weight_(weight), valuehat_(valuehat) {}
    
    SignalLikelihood(const SignalData& data) :
       value_(data.get_value()), weight_(data.get_weight()), valuehat_(data.get_valuehat()) {}
    
    NumericVector get_chi_square(double LB, double UB) const {
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> value_r = as<std::vector<double> >(value_);
        NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
        
        const NumericVector chisq = weight_ * SQUARE(valuehat_ - (soft + eCprime));
        return chisq;
    }

private:
    NumericVector value_, weight_, valuehat_;
};

#endif
