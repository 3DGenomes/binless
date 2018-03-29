#include "fast_dataframe.hpp"
#include "graph_helpers.hpp"

namespace binless {
namespace fast {

Rcpp::DataFrame get_as_dataframe(const FastData<Signal>& data, const BiasEstimator& bias, const DecayEstimator& dec, double tol_val) {
  //bias, decay, signal with decay and exposures, and log_background matrix (w/ offset)
  std::vector<double> biasmat,decaymat,signal,binless,background;
  biasmat.reserve(data.get_N());
  binless.reserve(data.get_N());
  std::vector<unsigned> dname(data.get_name()), nobs(data.get_nobs());
  std::vector<double> log_signal(data.get_log_signal()), exposures(data.get_exposures());
  Eigen::VectorXd log_biases = bias.get_data_estimate();
  Eigen::VectorXd log_decay = dec.get_data_estimate();
  for (unsigned i=0; i<data.get_N(); ++i) {
    double bipbj = log_biases(i);
    biasmat.push_back(std::exp(bipbj));
    double ldec = log_decay(i);
    decaymat.push_back(std::exp(ldec));
    double lsig = log_signal[i];
    signal.push_back(std::exp(lsig));
    unsigned name = dname[i]-1;
    double exposure = exposures[name];
    binless.push_back(std::exp(ldec + lsig));
    background.push_back(nobs[i]*std::exp(bipbj + ldec + exposure));
  }
  //build patchnos
  std::vector<double> patchnos;
  std::vector<unsigned> dbin1(data.get_bin1()), dbin2(data.get_bin2());
  patchnos.reserve(data.get_N());
  for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
    std::vector<unsigned> d_bin1(dbin1.cbegin() + dset*data.get_ncells(), dbin1.cbegin() + (dset+1)*data.get_ncells());
    std::vector<unsigned> d_bin2(dbin2.cbegin() + dset*data.get_ncells(), dbin2.cbegin() + (dset+1)*data.get_ncells());
    std::vector<double> d_beta(log_signal.cbegin() + dset*data.get_ncells(), log_signal.cbegin() + (dset+1)*data.get_ncells());
    IntegerVector patchno = get_patch_numbers(data.get_nbins(), tol_val, Rcpp::wrap(d_bin1), Rcpp::wrap(d_bin2), Rcpp::wrap(d_beta));
    for (const double v : patchno) patchnos.push_back(v);
  }
  //return
  return Rcpp::DataFrame::create(_["name"]=data.get_name_factor(),
                                 _["bin1"]=data.get_bin1_factor(),
                                 _["pos1"]=data.get_pos1(),
                                 _["bin2"]=data.get_bin2_factor(),
                                 _["pos2"]=data.get_pos2(),
                                 _["observed"]=data.get_observed(),
                                 _["nobs"]=data.get_nobs(),
                                 _["distance"]=data.get_distance(),
                                 _["background"]=background,
                                 _["biasmat"]=biasmat,
                                 _["decaymat"]=decaymat,
                                 _["signal"]=signal,
                                 _["binless"]=binless,
                                 _["phihat"]=data.get_signal_phihat(),
                                 _["weight"]=data.get_signal_weights(),
                                 _["patchno"]=patchnos);
}


Rcpp::DataFrame get_as_dataframe(const FastData<Difference>& data, double tol_val) {
  std::vector<double> log_difference = data.get_log_difference();
  std::vector<double> difference;
  difference.reserve(log_difference.size());
  for (const double v : log_difference) difference.push_back(std::exp(v));
  //build patchnos
  std::vector<unsigned> dbin1(data.get_bin1()), dbin2(data.get_bin2());
  std::vector<double> patchnos;
  patchnos.reserve(data.get_N());
  for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
    std::vector<unsigned> d_bin1(dbin1.cbegin() + dset*data.get_ncells(), dbin1.cbegin() + (dset+1)*data.get_ncells());
    std::vector<unsigned> d_bin2(dbin2.cbegin() + dset*data.get_ncells(), dbin2.cbegin() + (dset+1)*data.get_ncells());
    std::vector<double> d_beta(log_difference.cbegin() + dset*data.get_ncells(), log_difference.cbegin() + (dset+1)*data.get_ncells());
    IntegerVector patchno = get_patch_numbers(data.get_nbins(), tol_val, Rcpp::wrap(d_bin1), Rcpp::wrap(d_bin2), Rcpp::wrap(d_beta));
    for (const double v : patchno) patchnos.push_back(v);
  }
  return Rcpp::DataFrame::create(_["name"]=data.get_name_factor(),
                                 _["bin1"]=data.get_bin1_factor(),
                                 _["pos1"]=data.get_pos1(),
                                 _["bin2"]=data.get_bin2_factor(),
                                 _["pos2"]=data.get_pos2(),
                                 _["observed"]=data.get_observed(),
                                 _["nobs"]=data.get_nobs(),
                                 _["difference"]=difference,
                                 _["deltahat"]=data.get_deltahat(),
                                 _["weight"]=data.get_difference_weights(),
                                 _["phi_ref"]=data.get_phi_ref(),
                                 _["patchno"]=patchnos);
}

}
}

