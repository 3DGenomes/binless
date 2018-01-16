#include "fast_dataframe.hpp"


namespace binless {
namespace fast {

Rcpp::DataFrame get_as_dataframe(const FastData<Signal>& data, const DecayEstimator& dec) {
  //bias, decay, signal with decay and exposures, and log_background matrix (w/ offset)
  std::vector<double> biasmat,binless,log_background;
  biasmat.reserve(data.get_N());
  binless.reserve(data.get_N());
  std::vector<unsigned> dname(data.get_name()), dbin1(data.get_bin1()), dbin2(data.get_bin2());
  std::vector<double> log_biases(data.get_log_biases()), log_signal(data.get_log_signal()), exposures(data.get_exposures());
  Eigen::VectorXd decaymat = dec.get_data_estimate();
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
    unsigned bin2 = dbin2[i]-1;
    double bi = log_biases[bin1];
    double bj = log_biases[bin2];
    biasmat.push_back(bi + bj);
    double decay = decaymat(i);
    double signal = log_signal[i];
    unsigned name = dname[i]-1;
    double exposure = exposures[name];
    binless.push_back(decay + signal + exposure);
    log_background.push_back(bi + bj + decay + exposure);
  }
  return Rcpp::DataFrame::create(_["name"]=dname,
                                 _["bin1"]=dbin1,
                                 _["bin2"]=dbin2,
                                 _["observed"]=data.get_observed(),
                                 _["distance"]=data.get_distance(),
                                 _["log_background"]=log_background,
                                 _["log_biases"]=biasmat,
                                 _["log_decay"]=decaymat,
                                 _["log_signal"]=log_signal,
                                 _["log_binless"]=binless,
                                 _["phihat"]=data.get_signal_phihat(),
                                 _["weights"]=data.get_signal_weights() );
}


Rcpp::DataFrame get_as_dataframe(const FastData<Difference>& data) {
  return Rcpp::DataFrame::create(_["name"]=data.get_name(),
                                 _["bin1"]=data.get_bin1(),
                                 _["bin2"]=data.get_bin2(),
                                 _["observed"]=data.get_observed(),
                                 _["log_difference"]=data.get_log_difference(),
                                 _["deltahat"]=data.get_deltahat(),
                                 _["weights"]=data.get_difference_weights(),
                                 _["phi_ref"]=data.get_phi_ref());
}

}
}

