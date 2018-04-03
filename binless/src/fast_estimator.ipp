//here, each estimate is is log ( sum_i observed / sum_i expected ) with i representing a summary bin
template<class EstimatorImpl>
void Estimator<EstimatorImpl>::set_poisson_lsq_summary(const std::vector<double>& log_expected_std, const FastSignalData& data, double pseudocount) {
  //compute observed data
  auto observed_std = data.get_observed();
  auto nobs_std = data.get_nobs();
  Eigen::VectorXd observed = Eigen::VectorXd::Zero(observed_std.size());
  Eigen::VectorXd nobs = Eigen::VectorXd::Zero(nobs_std.size());
  for (unsigned i=0; i<observed.rows(); ++i) observed(i) = observed_std[i]; // cast to double
  for (unsigned i=0; i<observed.rows(); ++i) nobs(i) = nobs_std[i]; // cast to double
  Eigen::VectorXd sum_obs = summarize(observed);
  //compute expected data
  const Eigen::Map<const Eigen::VectorXd> log_expected(log_expected_std.data(),log_expected_std.size());
  Eigen::VectorXd sum_exp = summarize((log_expected.array().exp()*nobs.array()).matrix());
  //compute estimate
  Eigen::VectorXd estimate = (pseudocount + sum_obs.array()/sum_exp.array()).log().matrix();
  //center estimate
  double avg = EstimatorImpl::get_settings().get_nobs().dot(estimate)/EstimatorImpl::get_settings().get_nobs().sum();
  estimate = (estimate.array() - avg).matrix();
  //compute weight
  Eigen::VectorXd weight = (estimate.array().exp()*EstimatorImpl::get_settings().get_nobs().array()).matrix();
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance phihat weight sum_obs sum_exp\n";
  Rcpp::Rcout << (Eigen::MatrixXd(estimate.rows(),5) << EstimatorImpl::get_settings().get_log_distance().array().exp().matrix(),
                  estimate, weight, sum_obs, sum_exp).finished();*/
  //report back
  EstimatorImpl::get_summary().set_phihat(estimate);
  EstimatorImpl::get_summary().set_weight(weight);
  Rcpp::Rcout <<"set_poisson_lsq_summary for decay " << estimate.rows() << " " << weight.rows() << "\n";
  Rcpp::Rcout << EstimatorImpl::get_summary().get_phihat().rows() << "\n";
}

template<class EstimatorImpl>
void Estimator<EstimatorImpl>::update_summary(const ResidualsPair& z) {
  //compute weight
  const Eigen::Map<const Eigen::VectorXd> weights(z.weights.data(),z.weights.size());
  Eigen::VectorXd weight_sum = summarize(weights);
  //compute phihat
  const Eigen::Map<const Eigen::VectorXd> residuals(z.residuals.data(),z.residuals.size());
  Eigen::VectorXd phihat = summarize(residuals.array() * weights.array()).matrix();
  Eigen::VectorXd estimate = get_binned_estimate();
  phihat = (phihat.array() / weight_sum.array()).matrix() + estimate;
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance phihat weight\n";
  Rcpp::Rcout << (Eigen::MatrixXd(phihat.rows(),3) << EstimatorImpl::get_settings().get_log_distance().array().exp().matrix(), phihat,
                  weight_sum).finished();*/
  EstimatorImpl::get_summary().set_phihat(phihat);
  EstimatorImpl::get_summary().set_weight(weight_sum);
}
