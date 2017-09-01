#ifndef RAW_DATA_HPP
#define RAW_DATA_HPP

#include <Rcpp.h>
#include <vector>

class RawData {
public:
    RawData(int nbins, double dispersion, const Rcpp::DataFrame& cts, const Rcpp::List& outliers) :
    nbins_(nbins), dispersion_(dispersion), cts_(cts), outliers_(outliers) {}
    
    int get_nbins() const { return nbins_; }
    
    double get_dispersion() const { return dispersion_; }
    
    Rcpp::DataFrame get_cts() const { return cts_; }
    
    Rcpp::List get_outliers() const { return outliers_; }
    
    Rcpp::DataFrame get_mat() const { return mat_; }
    void set_mat(const Rcpp::DataFrame& mat) { mat_ = mat; }
    
    std::vector<double> get_beta() const { return beta_; }
    void set_beta(std::vector<double> beta) { beta_ = beta; }
    
private:
    const int nbins_;
    const double dispersion_;
    const Rcpp::DataFrame cts_;
    const Rcpp::List outliers_;
    
    Rcpp::DataFrame mat_;
    
    std::vector<double> beta_;
    
    
};

class SignalRawData : public RawData {
public:
    SignalRawData(int nbins, double dispersion, const Rcpp::DataFrame& cts, const Rcpp::List& outliers) :
     RawData(nbins, dispersion, cts, outliers) {}
};


class DifferenceRawData : public RawData {
public:
    DifferenceRawData(int nbins, double dispersion, const Rcpp::DataFrame& cts, const Rcpp::DataFrame& ref,
            const Rcpp::List& outliers) :
      RawData(nbins, dispersion, cts, outliers), ref_(ref) {}
    
    Rcpp::DataFrame get_ref() const { return ref_; }
    
    std::vector<double> get_phi_ref() const { return phi_ref_; }
    void set_phi_ref(const std::vector<double>& phi_ref) { phi_ref_ = phi_ref; }
    
private:
    const Rcpp::DataFrame ref_;
    
    std::vector<double> phi_ref_;
    
};



#endif

