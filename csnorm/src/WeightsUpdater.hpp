#ifndef WEIGHTS_UPDATER_HPP
#define WEIGHTS_UPDATER_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp" //compute_phi_ref
#include "cts_to_mat.hpp" //cts_to_diff_mat
#include "RawData.hpp"
#include "BinnedData.hpp"
#include "Traits.hpp"

//policy class that updates IRLS data and weights given counts and a new beta
template<typename Calculation> class WeightsUpdater;

//This class takes care of the case of signal estimation
template<> class WeightsUpdater<Signal> {
public:
    WeightsUpdater(const RawData<Signal>& raw, BinnedData<Signal>& binned) : raw_(raw), binned_(binned) {}
    
    std::vector<double> get_betahat() const { return Rcpp::as<std::vector<double> >(binned_.get_phihat()); }
    
    std::vector<double> get_weight() const { return Rcpp::as<std::vector<double> >(binned_.get_weight()); }
    
    void setUp() {}

    void update(const std::vector<double>& beta_phi) {
        Rcpp::NumericVector beta_phi_r = wrap(beta_phi);
        cts_to_signal_mat(raw_, 0, beta_phi_r, binned_); //offset is held at zero since we pass unthresholded beta_phi
    }

private:
    const RawData<Signal>& raw_;
    BinnedData<Signal>& binned_;
};

//This class takes care of the case of difference estimation
template<> class WeightsUpdater<Difference> {
public:
    WeightsUpdater(const RawData<Difference>& raw, BinnedData<Difference>& binned) : raw_(raw), binned_(binned) {}
    
    std::vector<double> get_betahat() const { return Rcpp::as<std::vector<double> >(binned_.get_deltahat()); }
    
    std::vector<double> get_weight() const { return Rcpp::as<std::vector<double> >(binned_.get_weight()); }
    
    void setUp(const std::vector<double>& phi_ref, const std::vector<double>& beta_delta) {
        //store first estimate of phi_ref and compute mat
        Rcpp::NumericVector phi_ref_r = wrap(phi_ref);
        Rcpp::NumericVector beta_delta_r = wrap(beta_delta);
        cts_to_diff_mat(raw_, phi_ref_r, beta_delta_r, binned_);
    }
    
    void update(const std::vector<double>& beta_delta) {
        //update phi_ref based on new beta
        //use beta because we don't threshold
        Rcpp::NumericVector beta_delta_r = wrap(beta_delta);
        Rcpp::NumericVector phi_ref_r = compute_phi_ref(binned_, beta_delta_r);
        //compute new matrix of weights
        cts_to_diff_mat(raw_, phi_ref_r, beta_delta_r, binned_);
    }
    
private:
    const RawData<Difference>& raw_;
    BinnedData<Difference>& binned_;
};

#endif

