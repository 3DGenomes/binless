#include "FastData.hpp"

namespace binless {
namespace fast {

Rcpp::DataFrame FastData<Signal>::get_as_dataframe() const {
    //bias, decay, signal with decay and exposures, and log_background matrix (w/ offset)
    std::vector<double> biasmat,decaymat,binless,log_background;
    biasmat.reserve(get_N());
    decaymat.reserve(get_N());
    binless.reserve(get_N());
    std::vector<unsigned> dname(get_name()), dbin1(get_bin1()), dbin2(get_bin2());
    std::vector<double> log_biases(get_log_biases()), log_decay(get_log_decay()),
    log_signal(get_log_signal()), exposures(get_exposures());
    for (unsigned i=0; i<get_N(); ++i) {
        unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = dbin2[i]-1;
        double bi = log_biases[bin1];
        double bj = log_biases[bin2];
        biasmat.push_back(bi + bj);
        double decay = log_decay[bin2-bin1];
        decaymat.push_back(decay);
        double signal = log_signal[i];
        unsigned name = dname[i]-1;
        double exposure = exposures[name];
        binless.push_back(decay + signal + exposure);
        log_background.push_back(bi + bj + decay + exposure);
    }
    return Rcpp::DataFrame::create(_["name"]=dname,
                                   _["bin1"]=dbin1,
                                   _["bin2"]=dbin2,
                                   _["observed"]=get_observed(),
                                   _["log_background"]=log_background,
                                   _["log_biases"]=biasmat,
                                   _["log_decay"]=decaymat,
                                   _["log_signal"]=log_signal,
                                   _["log_binless"]=binless,
                                   _["phihat"]=get_signal_phihat(),
                                   _["weights"]=get_signal_weights() );
}

std::vector<double> FastData<Difference>::get_log_signal() const {
    std::vector<double> delta = get_log_difference();
    std::vector<double> log_signal;
    log_signal.reserve(get_N());
    for (unsigned i=0; i<get_N(); ++i) {
        log_signal[i] = phi_ref_[i] + delta[i];
    }
    return log_signal;
}
void FastData<Difference>::set_log_signal(const std::vector<double>& log_signal) {
    //build phi_ref
    std::vector<double> phi_ref;
    phi_ref.reserve(get_N());
    for (unsigned dset=0; dset<get_ndatasets(); ++dset) {
        phi_ref.insert(phi_ref.end(), log_signal.cbegin() + ref_*get_ncells(), log_signal.cbegin()+(ref_+1)*get_ncells());
    }
    //build delta
    std::vector<double> delta;
    delta.reserve(get_N());
    for (unsigned i=0; i<get_N(); ++i) {
        delta.push_back(log_signal[i] - phi_ref[i]);
    }
    //assign
    set_log_difference(delta);
    std::swap(phi_ref,phi_ref_);
}

Rcpp::DataFrame FastData<Difference>::get_as_dataframe() const {
    return Rcpp::DataFrame::create(_["name"]=get_name(),
                                   _["bin1"]=get_bin1(),
                                   _["bin2"]=get_bin2(),
                                   _["observed"]=get_observed(),
                                   _["log_difference"]=get_log_difference(),
                                   _["deltahat"]=get_deltahat(),
                                   _["weights"]=get_difference_weights(),
                                   _["phi_ref"]=phi_ref_);
}

}
}

