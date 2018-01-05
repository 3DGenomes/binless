template<typename Derived>
std::vector<double> get_log_expected(const FastData<Derived>& data, const DecayEstimator& dec) {
    //
    std::vector<double> log_expected;
    log_expected.reserve(data.get_N());
    std::vector<unsigned> dbin1 = data.get_bin1();
    std::vector<unsigned> dbin2 = data.get_bin2();
    std::vector<double> dlog_biases = data.get_log_biases();
    std::vector<double> dsignal = data.get_log_signal();
    std::vector<unsigned> dname = data.get_name();
    std::vector<double> dexposures = data.get_exposures();
    Eigen::VectorXd dlog_decay = dec.get_data_log_decay();
    //
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = dbin2[i]-1;
        double bi = dlog_biases[bin1];
        double bj = dlog_biases[bin2];
        double log_decay = dlog_decay(i);
        double signal = dsignal[i];
        unsigned name = dname[i]-1;
        double exposure = dexposures[name];
        log_expected.push_back(bi + bj + log_decay + signal + exposure);
    }
    return log_expected;
}
