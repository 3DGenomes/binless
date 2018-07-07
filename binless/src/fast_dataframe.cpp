#include "fast_dataframe.hpp"
#include "graph_helpers.hpp"
#include "csvfile.hpp"
#include <algorithm>

namespace binless {
namespace fast {

Rcpp::DataFrame get_as_dataframe(const FastData<Signal>& data, const ExposureEstimator& expo, const BiasEstimator& bias,
                                 const DecayEstimator& dec, const NumericVector lam1, double tol_val, bool compute_patchnos,
                                 const std::string csv_out) {
  //bias, decay, signal with decay and exposures, and log_background matrix (w/ offset)
  std::vector<double> biasmat,decaymat,signal,binless,background;
  std::vector<double> log_signal(data.get_log_signal());
  std::vector<unsigned> nobs(data.get_nobs()), name(data.get_name());
  
  Eigen::VectorXd exposures = expo.get_data_estimate();
  Eigen::VectorXd log_biases = bias.get_data_estimate();
  Eigen::VectorXd log_decay = dec.get_data_estimate();
  const Eigen::ArrayXd lambda1_vals = (lam1.size() == data.get_ndatasets())
    ? Rcpp::as<Eigen::ArrayXd>(lam1)
      : Eigen::ArrayXd::Constant(data.get_ndatasets(),lam1[0]);
  if(!csv_out.empty()) {
    std::vector<double> vec_signal_phihat = data.get_signal_phihat();
    std::vector<double> vec_signal_weights = data.get_signal_weights();
    try
    {
      csvfile csv(csv_out);
      biasmat.push_back(0);
      decaymat.push_back(0);
      signal.push_back(0);
      binless.push_back(0);
      background.push_back(0);
      // Header
      csv << "phihat" << "weight";
      csv << "biasmat" << "decaymat" << "signal" << "binless" << "background" << endrow;
      for (unsigned i=0; i<data.get_N(); ++i) {
        double bipbj = log_biases(i);
        biasmat[0] = std::exp(bipbj);
        double ldec = log_decay(i);
        decaymat[0] = std::exp(ldec);
        double l1val = lambda1_vals(name[i]-1);
        double lsig = log_signal[i];
        lsig = (lsig > 0) ? std::max(0., lsig - l1val) : std::min(0., lsig + l1val);
        signal[0] = std::exp(lsig);
        double exposure = exposures(i);
        binless[0] = std::exp(ldec + lsig);
        background[0] = nobs[i]*std::exp(bipbj + ldec + exposure);
        csv << vec_signal_phihat[i] << vec_signal_weights[i];
        csv << biasmat[0] << decaymat[0] << signal[0] << binless[0] << background[0] << endrow;
      }
    }
    catch (const std::exception &ex)
    {
      Rcpp::Rcout << "Exception creating csv file: " << ex.what() << std::endl;
    }
    return 1;
  } else {
    biasmat.reserve(data.get_N());
    binless.reserve(data.get_N()); 
  
    for (unsigned i=0; i<data.get_N(); ++i) {
      double bipbj = log_biases(i);
      biasmat.push_back(std::exp(bipbj));
      double ldec = log_decay(i);
      decaymat.push_back(std::exp(ldec));
      double l1val = lambda1_vals(name[i]-1);
      double lsig = log_signal[i];
      lsig = (lsig > 0) ? std::max(0., lsig - l1val) : std::min(0., lsig + l1val);
      signal.push_back(std::exp(lsig));
      double exposure = exposures(i);
      binless.push_back(std::exp(ldec + lsig));
      background.push_back(nobs[i]*std::exp(bipbj + ldec + exposure));
    }
    //build patchnos
    std::vector<double> patchnos;
    patchnos.reserve(data.get_N());
    if (compute_patchnos) {
      std::vector<unsigned> dbin1(data.get_bin1()), dbin2(data.get_bin2());
      for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
        std::vector<unsigned> d_bin1(dbin1.cbegin() + dset*data.get_ncells(), dbin1.cbegin() + (dset+1)*data.get_ncells());
        std::vector<unsigned> d_bin2(dbin2.cbegin() + dset*data.get_ncells(), dbin2.cbegin() + (dset+1)*data.get_ncells());
        std::vector<double> d_beta(log_signal.cbegin() + dset*data.get_ncells(), log_signal.cbegin() + (dset+1)*data.get_ncells());
        IntegerVector patchno = get_patch_numbers(data.get_nbins(), tol_val, Rcpp::wrap(d_bin1), Rcpp::wrap(d_bin2), Rcpp::wrap(d_beta));
        for (const double v : patchno) patchnos.push_back(v);
      }
    } else {
      patchnos.assign(data.get_N(),-1.);
    }
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
  
  
}


Rcpp::DataFrame get_as_dataframe(const FastData<Difference>& data, const NumericVector lam1, double tol_val, bool compute_patchnos) {
  std::vector<double> log_difference = data.get_log_difference();
  std::vector<double> difference;
  difference.reserve(log_difference.size());
  for (const double v : log_difference) difference.push_back(std::exp(v));
  //build patchnos
  std::vector<double> patchnos;
  patchnos.reserve(data.get_N());
  if (compute_patchnos) {
    std::vector<unsigned> dbin1(data.get_bin1()), dbin2(data.get_bin2());
    for (unsigned dset=0; dset<data.get_ndatasets(); ++dset) {
      std::vector<unsigned> d_bin1(dbin1.cbegin() + dset*data.get_ncells(), dbin1.cbegin() + (dset+1)*data.get_ncells());
      std::vector<unsigned> d_bin2(dbin2.cbegin() + dset*data.get_ncells(), dbin2.cbegin() + (dset+1)*data.get_ncells());
      std::vector<double> d_beta(log_difference.cbegin() + dset*data.get_ncells(), log_difference.cbegin() + (dset+1)*data.get_ncells());
      IntegerVector patchno = get_patch_numbers(data.get_nbins(), tol_val, Rcpp::wrap(d_bin1), Rcpp::wrap(d_bin2), Rcpp::wrap(d_beta));
      for (const double v : patchno) patchnos.push_back(v);
    }
  } else {
    patchnos.assign(data.get_N(),-1.);
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

