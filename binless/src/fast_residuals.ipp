

template<typename FastData>
ResidualsPair get_residuals(const NormalDistribution& dist, const FastData& data, const ExposureEstimator& expo, const BiasEstimator& bias, const DecayEstimator& dec) {
    std::vector<double> residuals;
    std::vector<double> weights;
    residuals.reserve(data.get_N());
    weights.reserve(data.get_N());
    auto log_expected = get_log_expected(data, expo, bias, dec);
    auto observed = data.get_observed();
    auto nobs = data.get_nobs();
    for (unsigned i=0; i<data.get_N(); ++i) {
        if (observed[i]>0 && nobs[i]>0) {
            residuals.push_back( log(observed[i]/nobs[i]) - log_expected[i] );
            weights.push_back( nobs[i] );
        } else {
            residuals.push_back(  0 );
            weights.push_back( 0 );
        }
    }
    return ResidualsPair{residuals,weights};
}

template<typename FastData>
ResidualsPair get_residuals(const PoissonDistribution& dist, const FastData& data, const ExposureEstimator& expo, const BiasEstimator& bias, const DecayEstimator& dec) {
    std::vector<double> residuals;
    std::vector<double> weights;
    residuals.reserve(data.get_N());
    weights.reserve(data.get_N());
    auto log_expected = get_log_expected(data, expo, bias, dec);
    auto observed = data.get_observed();
    auto nobs = data.get_nobs();
    for (unsigned i=0; i<data.get_N(); ++i) {
        double expected_i = std::exp(log_expected[i]);
        residuals.push_back( observed[i]/(nobs[i] * expected_i) - 1 );
        weights.push_back( nobs[i] * expected_i );
    }
    return ResidualsPair{residuals,weights};
}

//residuals: negative binomial with log-link
template<typename FastData>
ResidualsPair get_residuals(const NegativeBinomialDistribution& dist, const FastData& data, const ExposureEstimator& expo, const BiasEstimator& bias, const DecayEstimator& dec) {
  std::vector<double> residuals;
  std::vector<double> weights;
  residuals.reserve(data.get_N());
  weights.reserve(data.get_N());
  auto log_expected = get_log_expected(data, expo, bias, dec);
  auto observed = data.get_observed();
  auto nobs = data.get_nobs();
  for (unsigned i=0; i<data.get_N(); ++i) {
    double expected_i = std::exp(log_expected[i]);
    residuals.push_back( observed[i]/(nobs[i] * expected_i) - 1);
    weights.push_back( nobs[i]/(1/expected_i + 1/dist.alpha) );
  }
  return ResidualsPair{residuals,weights};
}

