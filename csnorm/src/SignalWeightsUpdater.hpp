#ifndef SIGNAL_WEIGHTS_UPDATER_HPP
#define SIGNAL_WEIGHTS_UPDATER_HPP

#include <Rcpp.h>
#include <vector>

#include "cts_to_mat.hpp"
#include "RawData.hpp"

//policy class that updates IRLS data and weights given counts and a new beta
//This class takes care of the case of signal estimation
class SignalWeightsUpdater {
public:
    
    SignalWeightsUpdater(SignalRawData& data) : data_(data) {}
    
    void setUp() {} //for consistency with other WeightsUpdaters
    
    void update(const std::vector<double>& beta) {
        DataFrame mat = cts_to_signal_mat(data_.get_cts(), data_.get_nbins(), data_.get_dispersion(), beta, 0, data_.get_outliers());
        data_.set_beta(beta);
        data_.set_mat(mat);
    }
    
    std::vector<double> get_y() const { return Rcpp::as<std::vector<double> >(data_.get_mat()["phihat"]); }
    
    std::vector<double> get_w() const { return Rcpp::as<std::vector<double> >(data_.get_mat()["weight"]); }
    
    std::vector<int> get_bin1() const { return Rcpp::as<std::vector<int> > (data_.get_mat()["bin1"]); }
    
    std::vector<int> get_bin2() const { return Rcpp::as<std::vector<int> > (data_.get_mat()["bin2"]); }
    
    Rcpp::DataFrame get_mat() const { return data_.get_mat(); }
    
private:
    SignalRawData& data_;
};

#endif

