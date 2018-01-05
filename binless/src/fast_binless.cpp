#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_binless.hpp"
#include "fast_decay.hpp"
#include "fast_dataframe.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "spline.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {

//here, exposures are log ( sum_i observed / sum_i expected ) with i summed over datasets
std::vector<double> compute_poisson_lsq_exposures(const FastSignalData& data, const DecayEstimator& dec) {
  //get observed and expected
  std::vector<double> sum_obs(data.get_ndatasets(),0);
  std::vector<double> sum_exp(data.get_ndatasets(),0);
  auto log_expected = get_log_expected(data, dec);
  auto observed = data.get_observed();
  //sum them for each dataset
  auto names = data.get_name();
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = names[i]-1; //offset by 1 for vector indexing
    sum_obs[name] += observed[i];
    sum_exp[name] += std::exp(log_expected[i]);
  }
  //compute exposure
  std::vector<double> exposures;
  exposures.reserve(data.get_ndatasets());
  for (unsigned i=0; i<data.get_ndatasets(); ++i) {
    exposures.push_back(std::log(sum_obs[i]/sum_exp[i]));
  }
  return exposures;
}

//one IRLS iteration for exposures, with a poisson model
std::vector<double> step_exposures(const FastSignalData& data, const DecayEstimator& dec) {
  //get residuals
  ResidualsPair z = get_poisson_residuals(data, dec);
  //average by name
  std::vector<double> exposures(data.get_ndatasets(),0);
  std::vector<double> weightsums(data.get_ndatasets(),0);
  auto names = data.get_name();
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = names[i]-1; //offset by 1 for vector indexing
    exposures[name] += z.residuals[i]*z.weights[i];
    weightsums[name] += z.weights[i];
  }
  //add current estimates
  auto expo_ori = data.get_exposures();
  for (unsigned i=0; i<data.get_ndatasets(); ++i) {
    exposures[i] = exposures[i]/weightsums[i] + expo_ori[i];
  }
  return exposures;
}

//here, biases are log ( sum_i observed / sum_i expected ) with i summed over rows/columns
std::vector<double> compute_poisson_lsq_log_biases(const FastSignalData& data, const DecayEstimator& dec) {
  //get observed and expected data
  std::vector<double> sum_obs(data.get_nbins(),0);
  std::vector<double> sum_exp(data.get_nbins(),0);
  auto log_expected = get_log_expected(data, dec);
  auto observed = data.get_observed();
  //sum them along the rows/columns
  auto dbin1 = data.get_bin1();
  auto dbin2 = data.get_bin2();
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
    unsigned bin2 = dbin2[i]-1;
    double expected = std::exp(log_expected[i]);
    sum_obs[bin1] += observed[i];
    sum_exp[bin1] += expected;
    if (bin1!=bin2) {
      sum_obs[bin2] += observed[i];
      sum_exp[bin2] += expected;
    }
  }
  //compute log_bias
  std::vector<double> log_bias;
  log_bias.reserve(data.get_nbins());
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    log_bias.push_back(std::log(sum_obs[i]/sum_exp[i]));
  }
  //center log_bias
  double avg = std::accumulate(log_bias.begin(), log_bias.end(), 0.)/log_bias.size();
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    log_bias[i] -= avg;
  }
  return log_bias;
}

//one IRLS iteration for log biases, with a poisson model
std::vector<double> step_log_biases(const FastSignalData& data, const DecayEstimator& dec) {
  //get residuals
  ResidualsPair z = get_poisson_residuals(data, dec);
  //sum them along the rows/columns
  auto dbin1 = data.get_bin1();
  auto dbin2 = data.get_bin2();
  std::vector<double> log_biases(data.get_nbins(),0);
  std::vector<double> weightsums(data.get_nbins(),0);
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
    unsigned bin2 = dbin2[i]-1;
    double w = z.weights[i];
    double b = z.residuals[i]*w;
    log_biases[bin1] += b;
    weightsums[bin1] += w;
    if (bin1!=bin2) {
      log_biases[bin2] += b;
      weightsums[bin2] += w;
    }
  }
  //add current bias and compute weighted average
  double avg=0;
  double wsum=0;
  auto ori_log_biases = data.get_log_biases();
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    log_biases[i] = (weightsums[i]>0) ? log_biases[i]/weightsums[i] : 0;
    log_biases[i] += ori_log_biases[i];
    avg += log_biases[i]*weightsums[i];
    wsum += weightsums[i];
  }
  avg = avg / wsum;
  //no smoothing, subtract average and return
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    log_biases[i] -= avg;
  }
  return log_biases;
}

PrecisionPair get_precision(const std::vector<double>& weights, const std::vector<double>& weights_old) {
  double delta = std::abs(weights[0]-weights_old[0]);
  double maxval = weights[0];
  double minval = weights[0];
  const unsigned N = weights.size();
  for (unsigned i=1; i < N; ++i) {
    delta = std::max(std::abs(weights[i]-weights_old[i]), delta);
    minval = std::min(weights[i],minval);
    maxval = std::max(weights[i],maxval);
  }
  return PrecisionPair{delta,delta/(maxval-minval)};
}

std::vector<double> remove_signal_degeneracy(const FastSignalData& data) {
  auto dbin1 = data.get_bin1();
  auto dbin2 = data.get_bin2();
  auto dname = data.get_name();
  std::vector<double> log_signal = data.get_log_signal();
  double max_signal = *std::max_element(log_signal.begin(),log_signal.end());
  //get minimum signal per row and per counter diagonal
  std::vector<std::vector<double> > min_per_row(data.get_nbins(),std::vector<double>(data.get_ndatasets(),max_signal));
  std::vector<std::vector<double> > min_per_diag(data.get_nbins(),std::vector<double>(data.get_ndatasets(),max_signal));
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = dname[i]-1;
    unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
    unsigned bin2 = dbin2[i]-1;
    min_per_diag[bin2-bin1][name] = std::min(min_per_diag[bin2-bin1][name], log_signal[i]);
    min_per_row[bin1][name] = std::min(min_per_row[bin1][name], log_signal[i]);
    if (bin1!=bin2) {
      min_per_row[bin2][name] = std::min(min_per_row[bin2][name], log_signal[i]);
    }
  }
  //remove from each signal row and counter diagonal
  double max_adjust;
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = dname[i]-1;
    unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
    unsigned bin2 = dbin2[i]-1;
    double adjust = std::max(std::max(min_per_row[bin1][name],min_per_row[bin2][name]),min_per_diag[bin2-bin1][name]);
    log_signal[i] = std::max(log_signal[i]-adjust,0.);
    max_adjust = (i==0) ? adjust : std::max(adjust,max_adjust);
  }
  return log_signal;
}

std::vector<double> shift_signal(const FastSignalData& data) {
  auto dname = data.get_name();
  std::vector<double> log_signal = data.get_log_signal();
  double max_signal = *std::max_element(log_signal.begin(),log_signal.end());
  //get minimum signal per dataset
  std::vector<double> min_per_dset(data.get_ndatasets(),max_signal);
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = dname[i]-1;
    min_per_dset[name] = std::min(min_per_dset[name], log_signal[i]);
  }
  //shift signals
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = dname[i]-1;
    log_signal[i] = std::max(log_signal[i]-min_per_dset[name],0.);
  }
  return log_signal;
}

List binless(const DataFrame obs, unsigned nbins, double lam2, unsigned ngibbs, double tol_val, unsigned bg_steps) {
  //initialize return values, exposures and fused lasso optimizer
  Rcpp::Rcout << "init\n";
  FastSignalData out(obs, nbins);
  DecayConfig conf(tol_val);
  DecayEstimator dec(out, conf);
  //
  out.set_exposures(compute_poisson_lsq_exposures(out, dec));
  //
  dec.set_poisson_lsq_summary(out);
  dec.update_params();
  //
  out.set_log_biases(compute_poisson_lsq_log_biases(out, dec));
  double current_tol_val = 1.;
  std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(out.get_ndatasets(),
                                                             FusedLassoGaussianEstimator<GFLLibrary>(nbins, current_tol_val/20.));
  //std::vector<DataFrame> diagnostics;
  for (unsigned step=1; step <= ngibbs; ++step) {
    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "step " << step;
    std::vector<double> expected,old_expected;
    //compute biases
    auto biases = step_log_biases(out, dec);
    out.set_log_biases(biases);
    //compute decay
    if (step <= bg_steps) old_expected = get_log_expected(out, dec);
    dec.step_log_decay(out);
    if (step <= bg_steps) expected = get_log_expected(out, dec);
    //compute signal
    if (step > bg_steps) {
      old_expected = get_log_expected(out, dec);
      auto signal = step_signal(out, dec, flos, lam2);
      out.set_log_signal(signal.beta);
      out.set_signal_phihat(signal.phihat);
      out.set_signal_weights(signal.weights);
      expected = get_log_expected(out, dec);
    }
    //compute precision and convergence
    auto precision = get_precision(expected,old_expected);
    Rcpp::Rcout << " : reached precision abs = " << precision.abs << " rel = " << precision.rel << "\n";
      bool converged = precision.rel <= tol_val && current_tol_val <= tol_val*1.01;
      if (converged && step <= bg_steps) {
        converged = false;
        bg_steps = step; //force signal computation at next step
      }
      //update tolerance
      current_tol_val = std::min(current_tol_val,std::max(tol_val, precision.abs));
      for (unsigned i=0; i<out.get_ndatasets(); ++i) flos[i].set_tol(current_tol_val/20);
      //shift signal
      if (converged || step == ngibbs) {
        auto adjust = shift_signal(out);
        out.set_log_signal(adjust);
      } else {
        auto adjust = remove_signal_degeneracy(out);
        out.set_log_signal(adjust);
      }
      //compute exposures
      auto exposures = step_exposures(out, dec);
      out.set_exposures(exposures);
      //diagnostics.push_back(get_as_dataframe(out,dec));
      if (converged) {
        Rcpp::Rcout << "converged\n";
        break;
      }
  }
  Rcpp::Rcout << "done\n";
  //finalize and return
  return Rcpp::List::create(_["mat"]=get_as_dataframe(out,dec), _["log_biases"]=out.get_log_biases(),
                            _["beta_diag"]=dec.get_beta_diag(), _["exposures"]=out.get_exposures(),
                            //_["diagnostics"]=diagnostics,
                            _["nbins"]=nbins);
}

Rcpp::List binless_eval_cv(const List obs, const NumericVector lam2, unsigned group, double tol_val) {
  //read normalization data
  Rcpp::Rcout << "init\n";
  const unsigned nbins = obs["nbins"];
  const Rcpp::DataFrame mat = Rcpp::as<Rcpp::DataFrame>(obs["mat"]);
  FastSignalData out(mat, nbins);
  if (out.get_ndatasets() != 1) Rcpp::stop("Provide only 1 dataset!");
  auto signal_ori = Rcpp::as<std::vector<double> >(mat["signal"]);
  out.set_log_signal(signal_ori); //fills-in phi_ref and delta
  //
  DecayConfig conf(tol_val);
  DecayEstimator dec(out, conf);
  auto beta_diag = Rcpp::as<Eigen::VectorXd >(obs["beta_diag"]);
  dec.set_beta_diag(beta_diag);
  //
  auto log_biases = Rcpp::as<std::vector<double> >(obs["log_biases"]);
  out.set_log_biases(log_biases);
  auto exposures = Rcpp::as<std::vector<double> >(obs["exposures"]);
  out.set_exposures(exposures);
  const double converge = tol_val/20.;
  std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(1, FusedLassoGaussianEstimator<GFLLibrary>(nbins, converge));
  //compute cv
  Rcpp::Rcout << "compute\n";
  std::vector<DataFrame> diagnostics;
  for (unsigned i=0; i<lam2.size(); ++i) {
    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "lambda2= " << lam2[i] << "\n";
    out.set_log_signal(signal_ori); //fills-in phi_ref and delta
    auto signal = step_signal(out, dec, flos, lam2[i], group);
    out.set_log_signal(signal.beta);
    out.set_signal_phihat(signal.phihat);
    out.set_signal_weights(signal.weights);
    diagnostics.push_back(get_as_dataframe(out,dec));
  }
  //finalize and return
  Rcpp::Rcout << "done\n";
  return Rcpp::wrap(diagnostics);
}

Rcpp::DataFrame binless_difference(const List obs, double lam2, unsigned ref, double tol_val) {
  //read normalization data
  Rcpp::Rcout << "init\n";
  const unsigned nbins = obs["nbins"];
  const Rcpp::DataFrame mat = Rcpp::as<Rcpp::DataFrame>(obs["mat"]);
  FastDifferenceData out(mat, nbins, ref);
  auto signal = Rcpp::as<std::vector<double> >(mat["log_signal"]);
  out.set_log_signal(signal); //fills-in phi_ref and delta
  //
  DecayConfig conf(tol_val);
  DecayEstimator dec(out, conf);
  auto beta_diag = Rcpp::as<Eigen::VectorXd >(obs["beta_diag"]);
  dec.set_beta_diag(beta_diag);
  //
  auto log_biases = Rcpp::as<std::vector<double> >(obs["log_biases"]);
  out.set_log_biases(log_biases);
  auto exposures = Rcpp::as<std::vector<double> >(obs["exposures"]);
  out.set_exposures(exposures);
  const double converge = tol_val/20.;
  std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(out.get_ndatasets(),
                                                             FusedLassoGaussianEstimator<GFLLibrary>(nbins, converge));
  //compute differences
  Rcpp::Rcout << "compute\n";
  auto diff = step_difference(out, dec, flos, lam2, ref);
  out.set_log_difference(diff.delta);
  out.set_deltahat(diff.deltahat);
  out.set_difference_weights(diff.weights);
  out.set_phi_ref(diff.phi_ref);
  //finalize and return
  Rcpp::Rcout << "done\n";
  return get_as_dataframe(out);
}

}
}