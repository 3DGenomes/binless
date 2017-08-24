#ifndef SIGNAL_WEIGHTS_UPDATER_HPP
#define SIGNAL_WEIGHTS_UPDATER_HPP

#include <Rcpp.h>
#include <vector>

#include "cts_to_mat.hpp"

//policy class that updates IRLS data and weights given counts and a new beta
//This class takes care of the case of signal estimation
class SignalWeightsUpdater {
public:
    
    SignalWeightsUpdater(unsigned nrows, double dispersion,
                         const Rcpp::DataFrame& cts, const Rcpp::List& outliers) :
       nrows_(nrows), dispersion_(dispersion), cts_(cts), outliers_(outliers) {}
    
    void setUp() {} //for consistency with other WeightsUpdaters
    
    void update(const std::vector<double>& beta) {
        mat_ = cts_to_signal_mat(cts_, nrows_, dispersion_, beta, 0, outliers_);
    }
    
    std::vector<double> get_y() const { return Rcpp::as<std::vector<double> >(mat_["phihat"]); }
    
    std::vector<double> get_w() const { return Rcpp::as<std::vector<double> >(mat_["weight"]); }
    
    std::vector<int> get_bin1() const { return Rcpp::as<std::vector<int> > (mat_["bin1"]); }
    
    std::vector<int> get_bin2() const { return Rcpp::as<std::vector<int> > (mat_["bin2"]); }
    
    Rcpp::DataFrame get_mat() const { return mat_; }
    
private:
    const unsigned nrows_;
    const double dispersion_;
    const Rcpp::DataFrame cts_;
    const Rcpp::List outliers_;
    
    Rcpp::DataFrame mat_;
    
};

#endif

