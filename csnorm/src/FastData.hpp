#ifndef FAST_DATA_HPP
#define FAST_DATA_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Traits.hpp"

template<typename Derived>
class FastDataCore {
public:
    FastDataCore(const DataFrame& obs, unsigned nbins);
    
    unsigned get_N() const { return N_; }
    unsigned get_nbins() const { return nbins_; }
    unsigned get_ncells() const { return ncells_; }
    unsigned get_ndatasets() const { return ndatasets_; }
    std::vector<unsigned> get_name() const { return name_; }
    std::vector<unsigned> get_bin1() const { return bin1_; }
    std::vector<unsigned> get_bin2() const { return bin2_; }
    std::vector<unsigned> get_observed() const { return observed_; }
    
    std::vector<double> get_log_biases() const { return log_biases_; }
    void set_log_biases(const std::vector<double>& log_biases) { log_biases_ = log_biases; }
    
    std::vector<double> get_log_decay() const { return log_decay_; }
    void set_log_decay(const std::vector<double>& log_decay) { log_decay_ = log_decay; }
    
    std::vector<double> get_beta() const { return beta_; }
    void set_beta(const std::vector<double>& beta) { beta_ = beta; }
    
    std::vector<double> get_betahat() const { return betahat_; }
    void set_betahat(const std::vector<double>& betahat) { betahat_ = betahat; }
    
    std::vector<double> get_weights() const { return weights_; }
    void set_weights(const std::vector<double>& weights) { weights_ = weights; }
    
    std::vector<double> get_exposures() const { return exposures_; }
    void set_exposures(const std::vector<double>& exposures) { exposures_ = exposures; }
    
    std::vector<double> get_log_expected() const;
    
private:
    const std::vector<unsigned> name_,bin1_,bin2_,observed_; //N(N+1)/2
    const unsigned nbins_, ncells_, ndatasets_, N_;
    std::vector<double> log_biases_,log_decay_; //N
    std::vector<double> beta_,betahat_,weights_; //N(N+1)/2
    std::vector<double> exposures_; //ndatasets
};

template<typename> class FastData {};
typedef FastData<Signal> FastSignalData;
typedef FastData<Difference> FastDifferenceData;

template<> class FastData<Signal> : public FastDataCore<FastData<Signal> > {
public:
    FastData(const DataFrame& obs, unsigned nbins) : FastDataCore(obs,nbins) {}
    
    std::vector<double> get_log_signal() const { return get_beta(); }
    void set_log_signal(const std::vector<double>& log_signal) { set_beta(log_signal); }
    
    std::vector<double> get_signal_phihat() const { return get_betahat(); }
    void set_signal_phihat(const std::vector<double>& phihat) { set_betahat(phihat); }
    
    std::vector<double> get_signal_weights() const { return get_weights(); }
    void set_signal_weights(const std::vector<double>& weights) { set_weights(weights); }

    Rcpp::DataFrame get_as_dataframe() const;

};

template<> class FastData<Difference> : public FastDataCore<FastData<Difference> > {
public:
    FastData(const DataFrame& obs, unsigned nbins, unsigned ref) :
     FastDataCore(obs,nbins), phi_ref_(std::vector<double>(get_N(),0)), ref_(ref-1) {
        std::vector<unsigned> name = get_name();
        if (*std::max_element(name.begin(), name.end()) < ref || ref < 1) Rcpp::stop("Ref is not in the allowed range");
    }
    
    std::vector<double> get_log_difference() const { return get_beta(); }
    void set_log_difference(const std::vector<double>& delta) { set_beta(delta); }
    
    std::vector<double> get_deltahat() const { return get_betahat(); }
    void set_deltahat(const std::vector<double>& deltahat) { set_betahat(deltahat); }
    
    std::vector<double> get_difference_weights() const { return get_weights(); }
    void set_difference_weights(const std::vector<double>& weights) { set_weights(weights); }
    
    std::vector<double> get_phi_ref() const { return phi_ref_; }
    void set_phi_ref(const std::vector<double>& phi_ref) { phi_ref_ = phi_ref; }
    
    std::vector<double> get_log_signal() const;
    void set_log_signal(const std::vector<double>& log_signal);
    
    Rcpp::DataFrame get_as_dataframe() const;
    
private:
    std::vector<double> phi_ref_;
    unsigned ref_;
};

#include "FastData.ipp"

#endif
