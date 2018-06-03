
template<typename Lasso>
SignalTriplet step_signal(const FastSignalData& data, const ResidualsPair& z,
                          std::vector<Lasso>& flos, const Eigen::ArrayXd& lam2, unsigned group) {
    //build signal matrix
    auto dlog_signal = data.get_log_signal();
    std::vector<double> phihat,weights;
    std::vector<unsigned> nobs = data.get_nobs();
    phihat.reserve(data.get_N());
    weights.reserve(data.get_N());
    auto dbin1 = data.get_bin1();
    auto dbin2 = data.get_bin2();
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned cvgroup = 1 + ( (dbin2[i]+dbin1[i]) % 2 ); // 2 cv groups in checkerboard pattern
        if (cvgroup==group) {
            phihat.push_back(-100); //has influence if lam2==0
            weights.push_back(0);
        } else {
            phihat.push_back(z.residuals[i]+dlog_signal[i]);
            weights.push_back(z.weights[i]);
        }
    }
    //compute signal on each dataset
    std::vector<double> signal;
    signal.reserve(data.get_N());
    for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
        //run fused lasso
        std::vector<double> y(phihat.cbegin() + dset*data.get_ncells(), phihat.cbegin() + (dset+1)*data.get_ncells());
        std::vector<double> wt(weights.cbegin() + dset*data.get_ncells(), weights.cbegin() + (dset+1)*data.get_ncells());
        std::vector<unsigned> no(nobs.cbegin() + dset*data.get_ncells(), nobs.cbegin() + (dset+1)*data.get_ncells());
        flos[dset].optimize(y, wt, lam2[dset]);
        //compute weighted average
        std::vector<double> beta = flos[dset].get();
        double avg=0;
        double wsum=0;
        for (unsigned i=0; i<data.get_ncells(); ++i) {
          avg += beta[i]*no[i];
          wsum += no[i];
        }
        avg = avg / wsum;
        //subtract and store
        for (unsigned i=0; i<data.get_ncells(); ++i) {
            beta[i] = beta[i] - avg;
        }
        signal.insert(signal.end(), beta.begin(), beta.end());
    }
    return SignalTriplet{phihat,weights,signal};
}


template<typename Lasso>
DifferenceQuadruplet step_difference(const FastDifferenceData& data, const ResidualsPair& z,
                                     std::vector<Lasso>& flos, const Eigen::ArrayXd& lam2, unsigned ref) {
    //build difference matrices
    auto dsignal = data.get_log_signal();
    auto dphi_ref = data.get_phi_ref();
    std::vector<double> phihat,deltahat,weights;
    phihat.reserve(data.get_N());
    deltahat.reserve(data.get_N());
    weights.reserve(data.get_N());
    for (unsigned i=0; i<data.get_N(); ++i) {
        phihat.push_back(z.residuals[i]+dsignal[i]);
        deltahat.push_back(phihat.back() - dphi_ref[i]);
        weights.push_back(z.weights[i]);
    }
    std::vector<double> phihat_ref,weights_ref;
    phihat_ref.reserve(data.get_N());
    weights_ref.reserve(data.get_N());
    for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
        phihat_ref.insert(phihat_ref.end(), phihat.cbegin() + (ref-1)*data.get_ncells(), phihat.cbegin()+ref*data.get_ncells());
        weights_ref.insert(weights_ref.end(), weights.cbegin() + (ref-1)*data.get_ncells(), weights.cbegin()+ref*data.get_ncells());
    }
    //compute difference on each dataset
    std::vector<double> delta;
    delta.reserve(data.get_N());
    bool ref_seen = false;
    for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
        if (dset==ref-1) { //subtract 1 for R factor vs c++ indices
            delta.resize(delta.size()+data.get_ncells(), 0);
            ref_seen = true;
        } else {
            //run fused lasso
            std::vector<double> y(deltahat.cbegin() + dset*data.get_ncells(), deltahat.cbegin() + (dset+1)*data.get_ncells());
            std::vector<double> wt(weights.cbegin() + dset*data.get_ncells(), weights.cbegin() + (dset+1)*data.get_ncells());
            unsigned idx = dset - (ref_seen ? 1 : 0);
            flos[idx].optimize(y, wt, lam2[idx]);
            std::vector<double> ldelta = flos[idx].get();
            delta.insert(delta.end(), ldelta.begin(), ldelta.end());
        }
    }
    //compute phi_ref
    std::vector<double> phi_ref;
    phi_ref.reserve(data.get_N());
    for (unsigned i=0; i<data.get_N(); ++i) {
        double val;
        if (weights_ref[i]==0 && weights[i]==0) {
            val = (phihat_ref[i]+phihat[i])/2;
        } else {
            val = (phihat_ref[i]*weights_ref[i] + (phihat[i]-delta[i])*weights[i])/(weights_ref[i]+weights[i]);
        }
        val = std::max(val,0.);
        phi_ref.push_back(val);
    }
    return DifferenceQuadruplet{deltahat,weights,delta,phi_ref};
}
