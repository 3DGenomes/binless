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

    void set_beta(const std::vector<double>& beta) { set_beta_phi(beta); }
    
    void update(const std::vector<double>& beta_phi) {
        set_beta_phi(beta_phi);
        cts_to_signal_mat(raw_, binned_); //offset is held at zero since we pass unthresholded beta_phi
    }

private:
    
    void set_beta_phi(const std::vector<double>& beta_phi) {
        binned_.set_beta_phi(wrap(beta_phi));
    }
   
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
        binned_.set_beta_delta(wrap(beta_delta));
        binned_.set_phi_ref(wrap(phi_ref));
        cts_to_diff_mat(raw_, binned_);
    }
   
    void set_beta(const std::vector<double>& beta) { set_beta_delta(beta); }
        
    void update(const std::vector<double>& beta_delta) {
        set_beta_delta(beta_delta);
        //compute new matrix of weights
        cts_to_diff_mat(raw_, binned_);
    }
    
private:
    
    void set_beta_delta(const std::vector<double>& beta_delta) {
        //update phi_ref based on new beta
        //use beta because we don't threshold
        binned_.set_beta_delta(wrap(beta_delta));
        binned_.set_phi_ref(compute_phi_ref(binned_, binned_.get_beta_delta()));
    }
    
    const RawData<Difference>& raw_;
    BinnedData<Difference>& binned_;
};

#endif

