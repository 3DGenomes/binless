#ifndef DIFFERENCE_WEIGHTS_UPDATER_HPP
#define DIFFERENCE_WEIGHTS_UPDATER_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp" //compute_phi_ref
#include "cts_to_mat.hpp" //cts_to_diff_mat
#include "RawData.hpp"
#include "BinnedData.hpp"

//policy class that updates IRLS data and weights given counts and a new beta
//This class takes care of the case of difference estimation
class DifferenceWeightsUpdater {
public:
    
    DifferenceWeightsUpdater(const DifferenceRawData& raw, BinnedData<Difference>& binned) : raw_(raw), binned_(binned) {}
    
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
    
    std::vector<double> get_betahat() const {
        std::vector<double> phihat = Rcpp::as<std::vector<double> >(binned_.get_phihat());
        const unsigned N = phihat.size();
        std::vector<double> y_r;
        y_r.reserve(N);
        std::vector<double> phi_ref = Rcpp::as<std::vector<double> >(binned_.get_phi_ref());
        for (int i=0; i<N; ++i) {
            y_r.push_back(phihat[i]-phi_ref[i]);
        }
        return y_r;
    }
    
    std::vector<double> get_weight() const {
        std::vector<double> w = Rcpp::as<std::vector<double> >(binned_.get_weight());
        return w;
    }
    
    std::vector<int> get_bin1() const { return Rcpp::as<std::vector<int> > (binned_.get_bin1()); }
    
    std::vector<int> get_bin2() const { return Rcpp::as<std::vector<int> > (binned_.get_bin2()); }
    
    std::vector<double> get_phi_ref() const { return Rcpp::as<std::vector<double> >(binned_.get_phi_ref()); }
    
private:
    const DifferenceRawData& raw_;
    BinnedData<Difference>& binned_;
};

#endif

