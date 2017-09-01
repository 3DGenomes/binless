#ifndef BINNED_DATA_HPP
#define BINNED_DATA_HPP

#include <Rcpp.h>

class BinnedData {
public:
    BinnedData(const Rcpp::NumericVector& value, const Rcpp::NumericVector& weight,
         const Rcpp::NumericVector& valuehat, const Rcpp::NumericVector& ncounts,
         const Rcpp::IntegerVector& patchno) :
    value_(value), weight_(weight), valuehat_(valuehat), ncounts_(ncounts), patchno_(patchno) {}
    
    Rcpp::NumericVector get_value() const { return value_; }
    
    Rcpp::NumericVector get_weight() const { return weight_; }
    
    Rcpp::NumericVector get_valuehat() const { return valuehat_; }
    
    Rcpp::IntegerVector get_patchno() const { return patchno_; }
    
private:
    const Rcpp::NumericVector value_, weight_, valuehat_, ncounts_;
    const Rcpp::IntegerVector patchno_;
};

class SignalBinnedData : public BinnedData {
public:
  SignalBinnedData(const Rcpp::NumericVector& value, const Rcpp::NumericVector& weight,
       const Rcpp::NumericVector& valuehat, const Rcpp::NumericVector& ncounts,
       const Rcpp::IntegerVector& patchno) :
   BinnedData(value, weight, valuehat, ncounts, patchno) {}
};

class DifferenceBinnedData : public BinnedData {
public:
    DifferenceBinnedData(const Rcpp::NumericVector& value, const Rcpp::NumericVector& weight,
                   const Rcpp::NumericVector& valuehat, const Rcpp::NumericVector& weight_ref,
                   const Rcpp::NumericVector& valuehat_ref, const Rcpp::NumericVector& ncounts,
                   const Rcpp::IntegerVector& patchno) :
      BinnedData(value, weight, valuehat, ncounts, patchno), weight_ref_(weight_ref), valuehat_ref_(valuehat_ref) {}
    
    Rcpp::NumericVector get_weight_ref() const { return weight_ref_; }
    
    Rcpp::NumericVector get_valuehat_ref() const { return valuehat_ref_; }
    
private:
    const Rcpp::NumericVector weight_ref_, valuehat_ref_;
};



#endif

