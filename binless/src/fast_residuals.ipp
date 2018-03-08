

//residuals: normal with log-link, 0 drops data
template<typename FastData>
ResidualsPair get_normal_residuals(const FastData& data, const DecayEstimator& dec) {
    std::vector<double> residuals;
    std::vector<double> weights;
    residuals.reserve(data.get_N());
    weights.reserve(data.get_N());
    auto log_expected = get_log_expected(data, dec);
    auto observed = data.get_observed();
    auto nobs = data.get_nobs();
    for (unsigned i=0; i<data.get_N(); ++i) {
        if (observed[i]>0) {
            residuals.push_back( log(observed[i]) - log_expected[i] );
            weights.push_back( nobs[i] );
        } else {
            residuals.push_back(  0 );
            weights.push_back( 0 );
        }
    }
    return ResidualsPair{residuals,weights};
}

//residuals: poisson with log-link
template<typename FastData>
ResidualsPair get_poisson_residuals(const FastData& data, const DecayEstimator& dec) {
    std::vector<double> residuals;
    std::vector<double> weights;
    residuals.reserve(data.get_N());
    weights.reserve(data.get_N());
    auto log_expected = get_log_expected(data, dec);
    auto observed = data.get_observed();
    auto nobs = data.get_nobs();
    for (unsigned i=0; i<data.get_N(); ++i) {
        double expected_i = std::exp(log_expected[i]);
        residuals.push_back( (observed[i]/expected_i) - 1);
        weights.push_back( expected_i*nobs[i] );
    }
    return ResidualsPair{residuals,weights};
}
