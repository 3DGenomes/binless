template<typename Derived>
std::vector<double> get_log_expected(const FastData<Derived>& data, const BiasEstimator& bias, const DecayEstimator& dec) {
    //
    std::vector<double> log_expected;
    log_expected.reserve(data.get_N());
    Eigen::VectorXd dlog_biases = bias.get_data_estimate();
    std::vector<double> dsignal = data.get_log_signal();
    std::vector<unsigned> dname = data.get_name();
    std::vector<double> dexposures = data.get_exposures();
    Eigen::VectorXd dlog_decay = dec.get_data_estimate();
    //
    for (unsigned i=0; i<data.get_N(); ++i) {
        double log_biases = dlog_biases(i);
        double log_decay = dlog_decay(i);
        double signal = dsignal[i];
        unsigned name = dname[i]-1;
        double exposure = dexposures[name];
        log_expected.push_back(log_biases + log_decay + signal + exposure);
    }
    return log_expected;
}
