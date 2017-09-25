#ifndef RAW_DATA_HPP
#define RAW_DATA_HPP

#include <Rcpp.h>
#include <vector>

#include "Traits.hpp"

class RawDataCore {
public:
    RawDataCore(int nbins, double dispersion, const Rcpp::DataFrame& cts, const Rcpp::List& outliers) :
    nbins_(nbins), dispersion_(dispersion), cts_(cts), outliers_(outliers) {}
    
    int get_nbins() const { return nbins_; }
    
    double get_dispersion() const { return dispersion_; }
    
    Rcpp::DataFrame get_cts() const { return cts_; }
    
    Rcpp::List get_outliers() const { return outliers_; }
    
private:
    const int nbins_;
    const double dispersion_;
    const Rcpp::DataFrame cts_;
    const Rcpp::List outliers_;
};

template<typename> class RawData {};

template<> class RawData<Signal> : public RawDataCore {
public:
    RawData(int nbins, double dispersion, const Rcpp::DataFrame& cts, const Rcpp::List& outliers) :
     RawDataCore(nbins, dispersion, cts, outliers) {}
};


template<> class RawData<Difference> : public RawDataCore {
public:
    RawData(int nbins, double dispersion, const Rcpp::DataFrame& cts, const Rcpp::DataFrame& ref,
            const Rcpp::List& outliers) :
      RawDataCore(nbins, dispersion, cts, outliers), ref_(ref) {}
    
    Rcpp::DataFrame get_ref() const { return ref_; }
    
private:
    const Rcpp::DataFrame ref_;
};



#endif

