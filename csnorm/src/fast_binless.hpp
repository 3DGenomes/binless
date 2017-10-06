#ifndef FAST_BINLESS_HPP
#define FAST_BINLESS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"

class FastData {
public:
    FastData(const DataFrame& obs, unsigned nbins);
    
    unsigned get_N() const { return N_; }
    unsigned get_nbins() const { return nbins_; }
    unsigned get_ncells() const { return ncells_; }
    unsigned get_ndatasets() const { return ndatasets_; }
    std::vector<unsigned> get_name() const { return name_; }
    std::vector<unsigned> get_bin1() const { return bin1_; }
    std::vector<unsigned> get_bin2() const { return bin2_; }
    std::vector<unsigned> get_observed() const { return observed_; }
     
    std::vector<double> get_log_biases() const { return log_biases_; }
    void set_log_biases(const std::vector<double> log_biases) { log_biases_ = log_biases; }
     
    std::vector<double> get_log_decay() const { return log_decay_; }
    void set_log_decay(const std::vector<double> log_decay) { log_decay_ = log_decay; }

    std::vector<double> get_log_signal() const { return log_signal_; }
    void set_log_signal(const std::vector<double> log_signal) { log_signal_ = log_signal; }
    
    std::vector<double> get_signal_phihat() const { return phihat_; }
    void set_signal_phihat(const std::vector<double> phihat) { phihat_ = phihat; }
    
    std::vector<double> get_signal_weights() const { return weights_; }
    void set_signal_weights(const std::vector<double> weights) { weights_ = weights; }

    std::vector<double> get_exposures() const { return exposures_; }
    void set_exposures(const std::vector<double> exposures) { exposures_ = exposures; }
    
    std::vector<double> get_log_expected() const;
    
    Rcpp::DataFrame get_as_dataframe() const;

private:
    const std::vector<unsigned> name_,bin1_,bin2_,observed_; //N(N+1)/2
    const unsigned nbins_, ncells_, ndatasets_, N_;
    std::vector<double> log_biases_,log_decay_; //N
    std::vector<double> log_signal_,phihat_,weights_; //N(N+1)/2
    std::vector<double> exposures_; //ndatasets
};


struct ResidualsPair { std::vector<double> residuals,weights; };
struct SignalTriplet { std::vector<double> phihat, weights, beta; };

ResidualsPair get_normal_residuals(const FastData& data);
ResidualsPair get_poisson_residuals(const FastData& data);
std::vector<double> fast_compute_exposures(const FastData& data);
std::vector<double> fast_compute_log_biases(const FastData& data);
std::vector<double> fast_compute_log_decay(const FastData& data);
template<typename Lasso>
SignalTriplet fast_compute_signal(const FastData& data, Lasso& flo, double lam2);

List fast_binless(const DataFrame obs, unsigned nbins, unsigned nouter, double lam2, double tol_val);

#include "fast_binless.ipp"

#endif

