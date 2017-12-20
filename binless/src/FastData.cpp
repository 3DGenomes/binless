#include "FastData.hpp"

namespace binless {
namespace fast {

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

}
}

