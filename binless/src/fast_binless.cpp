#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_binless.hpp"
#include "fast_decay.hpp"
#include "fast_bias_mean.hpp"
#include "fast_dataframe.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "spline.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {

//here, exposures are log ( sum_i observed / sum_i expected ) with i summed over datasets
std::vector<double> compute_poisson_lsq_exposures(const FastSignalData& data, const BiasEstimator& bias, const DecayEstimator& dec, double pseudocount) {
  //get observed and expected
  std::vector<double> sum_obs(data.get_ndatasets(),0);
  std::vector<double> sum_exp(data.get_ndatasets(),0);
  auto log_expected = get_log_expected(data, bias, dec);
  auto observed = data.get_observed();
  auto nobs = data.get_nobs();
  //sum them for each dataset
  auto names = data.get_name();
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned name = names[i]-1; //offset by 1 for vector indexing
    sum_obs[name] += observed[i];
    sum_exp[name] += std::exp(log_expected[i])*nobs[i];
  }
  //compute exposure
  std::vector<double> exposures;
  exposures.reserve(data.get_ndatasets());
  for (unsigned i=0; i<data.get_ndatasets(); ++i) {
    exposures.push_back(std::log(pseudocount + sum_obs[i]/sum_exp[i]));
  }
  return exposures;
}

//one IRLS iteration for exposures, with a poisson model
std::vector<double> step_exposures(const FastSignalData& data, const ResidualsPair& z) {
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

List binless(const DataFrame obs, unsigned nbins, double lam2, double alpha, unsigned ngibbs, double tol_val, unsigned bg_steps, unsigned free_decay) {
  //initialize return values, exposures and fused lasso optimizer
  Rcpp::Rcout << "init\n";
  NegativeBinomialDistribution nb_dist;
  Sampler<NegativeBinomialDistribution> nb_sampler(nb_dist);
  nb_sampler.init(alpha);
  FastSignalData out(obs, nbins);
  //
  DecayConfig dconf(tol_val, free_decay);
  DecayEstimator dec(out, dconf);
  //unsigned constraint_every = 0;
  //BiasConfig bconf(tol_val, constraint_every);
  BiasConfig bconf(nbins);
  BiasEstimator bias(out, bconf);
  //
  out.set_exposures(compute_poisson_lsq_exposures(out, bias, dec));
  //
  std::vector<double> expected;
  expected = get_log_expected(out, bias, dec);
  dec.set_poisson_lsq_summary(expected, out);
  dec.update_params();
  //
  expected = get_log_expected(out, bias, dec);
  bias.set_poisson_lsq_summary(expected, out);
  bias.update_params();
  double current_tol_val = 1.;
  std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(out.get_ndatasets(),
                                                             FusedLassoGaussianEstimator<GFLLibrary>(nbins, current_tol_val/20.));
  //std::vector<DataFrame> diagnostics;
  ResidualsPair z;
  for (unsigned step=1; step <= ngibbs; ++step) {
    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "step " << step;
    std::vector<double> old_expected;
    //compute biases
    z = get_residuals(nb_dist, out, bias, dec);
    bias.step_irls(z);
    //compute decay
    if (step <= bg_steps) old_expected = get_log_expected(out, bias, dec);
    z = get_residuals(nb_dist, out, bias, dec);
    dec.step_irls(z);
    if (step <= bg_steps) expected = get_log_expected(out, bias, dec);
    //compute signal
    if (step > bg_steps) {
      old_expected = get_log_expected(out, bias, dec);
      z = get_residuals(nb_dist, out, bias, dec);
      auto signal = step_signal(out, z, flos, lam2);
      out.set_log_signal(signal.beta);
      out.set_signal_phihat(signal.phihat);
      out.set_signal_weights(signal.weights);
      expected = get_log_expected(out, bias, dec);
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
      z = get_residuals(nb_dist, out, bias, dec);
      auto exposures = step_exposures(out, z);
      out.set_exposures(exposures);
      //diagnostics.push_back(get_as_dataframe(out,dec,tol_val));
      if (converged) {
        Rcpp::Rcout << "converged\n";
        break;
      }
  }
  Rcpp::Rcout << "done\n";
  //finalize and return
  return Rcpp::List::create(_["mat"]=get_as_dataframe(out, bias, dec, tol_val), _["log_biases"]=bias.get_binned_estimate(),
                            _["beta_diag"]=dec.get_beta(), _["exposures"]=out.get_exposures(),
                            _["log_signal"]=out.get_log_signal(),
                            //_["diagnostics"]=diagnostics,
                            _["nbins"]=nbins);
}

Rcpp::List binless_eval_cv(const List obs, const NumericVector lam2, double alpha, unsigned group, double tol_val) {
  //setup distribution
  NegativeBinomialDistribution nb_dist;
  Sampler<NegativeBinomialDistribution> nb_sampler(nb_dist);
  nb_sampler.init(alpha);
  //
  //read normalization data
  Rcpp::Rcout << "init\n";
  const unsigned nbins = obs["nbins"];
  const Rcpp::DataFrame mat = Rcpp::as<Rcpp::DataFrame>(obs["mat"]);
  FastSignalData out(mat, nbins);
  if (out.get_ndatasets() != 1) Rcpp::stop("Provide only 1 dataset!");
  auto signal_ori = Rcpp::as<std::vector<double> >(mat["signal"]);
  out.set_log_signal(signal_ori); //fills-in phi_ref and delta
  //
  DecayConfig conf(tol_val, 10000); //no need to pass free_diag as parameter because it is not used anyway
  DecayEstimator dec(out, conf);
  auto beta_diag = Rcpp::as<Eigen::VectorXd >(obs["beta_diag"]);
  dec.set_beta(beta_diag);
  unsigned constraint_every = 0;
  BiasConfig bconf(nbins); //bconf(tol_val, constraint_every);
  BiasEstimator bias(out, bconf);
  //
  bias.set_estimate(Rcpp::as<Eigen::VectorXd>(obs["log_biases"]));
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
    ResidualsPair z = get_residuals(nb_dist, out, bias, dec);
    auto signal = step_signal(out, z, flos, lam2[i], group);
    out.set_log_signal(signal.beta);
    out.set_signal_phihat(signal.phihat);
    out.set_signal_weights(signal.weights);
    diagnostics.push_back(get_as_dataframe(out, bias, dec, tol_val));
  }
  //finalize and return
  Rcpp::Rcout << "done\n";
  return Rcpp::wrap(diagnostics);
}

Rcpp::DataFrame binless_difference(const List obs, double lam2, unsigned ref, double alpha, double tol_val) {
  //setup distribution
  NegativeBinomialDistribution nb_dist;
  Sampler<NegativeBinomialDistribution> nb_sampler(nb_dist);
  nb_sampler.init(alpha);
  //
  //read normalization data
  Rcpp::Rcout << "init\n";
  const unsigned nbins = obs["nbins"];
  const Rcpp::DataFrame mat = Rcpp::as<Rcpp::DataFrame>(obs["mat"]);
  FastDifferenceData out(mat, nbins, ref);
  auto log_signal = Rcpp::as<std::vector<double> >(obs["log_signal"]);
  out.set_log_signal(log_signal); //fills-in phi_ref and delta
  //
  DecayConfig conf(tol_val, 10000); //no need to pass free_diag as parameter because it is not used anyway
  DecayEstimator dec(out, conf);
  auto beta_diag = Rcpp::as<Eigen::VectorXd >(obs["beta_diag"]);
  dec.set_beta(beta_diag);
  unsigned constraint_every = 0;
  BiasConfig bconf(nbins); //bconf(tol_val, constraint_every);
  BiasEstimator bias(out, bconf);
  //
  bias.set_estimate(Rcpp::as<Eigen::VectorXd>(obs["log_biases"]));
  auto exposures = Rcpp::as<std::vector<double> >(obs["exposures"]);
  out.set_exposures(exposures);
  const double converge = tol_val/20.;
  std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(out.get_ndatasets(),
                                                             FusedLassoGaussianEstimator<GFLLibrary>(nbins, converge));
  //compute differences
  Rcpp::Rcout << "compute\n";
  ResidualsPair z = get_residuals(nb_dist, out, bias, dec);
  auto diff = step_difference(out, z, flos, lam2, ref);
  out.set_log_difference(diff.delta);
  out.set_deltahat(diff.deltahat);
  out.set_difference_weights(diff.weights);
  out.set_phi_ref(diff.phi_ref);
  //finalize and return
  Rcpp::Rcout << "done\n";
  return get_as_dataframe(out,tol_val);
}

}
}