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
    
    SignalWeightsUpdater(const SignalRawData& raw, SignalBinnedData& binned) : raw_(raw), binned_(binned) {}
    
    void setUp() {} //for consistency with other WeightsUpdaters
    
    void update(const std::vector<double>& beta_phi) {
        Rcpp::NumericVector beta_phi_r = wrap(beta_phi);
        cts_to_signal_mat(raw_, 0, beta_phi_r, binned_); //offset is held at zero since we pass unthresholded beta_phi
    }
    
    std::vector<double> get_betahat() const { return Rcpp::as<std::vector<double> >(binned_.get_phihat()); }
    
    std::vector<double> get_weight() const { return Rcpp::as<std::vector<double> >(binned_.get_weight()); }
    
    std::vector<int> get_bin1() const { return Rcpp::as<std::vector<int> > (binned_.get_bin1()); }
    
    std::vector<int> get_bin2() const { return Rcpp::as<std::vector<int> > (binned_.get_bin2()); }
    
private:
    const SignalRawData& raw_;
    SignalBinnedData& binned_;
};

#endif

