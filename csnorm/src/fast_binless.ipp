
template<typename Lasso>
SignalTriplet fast_compute_signal(const FastData& data, Lasso& flo, double lam2) {
    //get residuals
    ResidualsPair z = get_poisson_residuals(data);
    //build signal matrix
    auto dlog_signal = data.get_log_signal();
    std::vector<double> phihat,weights;
    phihat.reserve(data.get_N());
    weights.reserve(data.get_N());
    //double wavg = double(std::accumulate(z.weights.begin(),z.weights.end(),0.))/z.weights.size();
    for (unsigned i=0; i<data.get_N(); ++i) {
        phihat.push_back(z.residuals[i]+dlog_signal[i]);
        weights.push_back(z.weights[i]); ///wavg);
    }
    //compute signal on each dataset
    std::vector<double> signal;
    signal.reserve(data.get_N());
    for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
        //run fused lasso
        std::vector<double> y(phihat.cbegin() + dset*data.get_ncells(), phihat.cbegin() + (dset+1)*data.get_ncells());
        std::vector<double> wt(weights.cbegin() + dset*data.get_ncells(), weights.cbegin() + (dset+1)*data.get_ncells());
        flo.optimize(y, wt, lam2);
        //subtract average and return
        std::vector<double> beta = flo.get();
        double avg = std::accumulate(beta.begin(),beta.end(),0.)/beta.size();
        for (unsigned i=0; i<data.get_ncells(); ++i) {
            beta[i] = beta[i] - avg;
        }
        signal.insert(signal.end(), beta.begin(), beta.end());
    }
    return SignalTriplet{phihat,weights,signal};
}
