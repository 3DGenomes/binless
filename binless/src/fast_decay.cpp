#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fast_decay.hpp"
#include "fast_residuals.hpp"
#include "spline.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {


//here, log decay is log ( sum_i observed / sum_i expected ) with i summed over counter diagonals
std::vector<double> compute_poisson_lsq_log_decay(const FastSignalData& data) {
  //get observed and expected data
  std::vector<double> sum_obs(data.get_nbins(),0);
  std::vector<double> sum_exp(data.get_nbins(),0);
  auto log_expected = data.get_log_expected();
  auto observed = data.get_observed();
  //sum them along the counter diagonals
  auto dbin1 = data.get_bin1();
  auto dbin2 = data.get_bin2();
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned dist = dbin2[i]-dbin1[i];
    double expected = std::exp(log_expected[i]);
    sum_obs[dist] += observed[i];
    sum_exp[dist] += expected;
  }
  //compute log_decay
  std::vector<double> log_decay;
  log_decay.reserve(data.get_nbins());
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    if (sum_obs[i]==0) {
      Rcpp::Rcout << "counter diag " << i << " is zero!\n";
      Rcpp::stop(" Aborting...");
    }
    log_decay.push_back(std::log(sum_obs[i]/sum_exp[i]));
  }
  //center log_decay
  double avg = std::accumulate(log_decay.begin(), log_decay.end(), 0.)/log_decay.size();
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    log_decay[i] -= avg;
  }
  return log_decay;
}

DecaySummary get_decay_summary(const FastSignalData& data) {
  //get residuals
  ResidualsPair z = get_poisson_residuals(data);
  //sum them along the diagonals
  auto dbin1 = data.get_bin1();
  auto dbin2 = data.get_bin2();
  std::vector<double> kappahat(data.get_nbins(),0);
  std::vector<double> weight(data.get_nbins(),0);
  for (unsigned i=0; i<data.get_N(); ++i) {
    unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
    unsigned bin2 = dbin2[i]-1;
    kappahat[bin2-bin1] += z.residuals[i]*z.weights[i];
    weight[bin2-bin1] += z.weights[i];
  }
  //add current bias and normalize
  auto log_decay = data.get_log_decay();
  std::vector<double> distance;
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    kappahat[i] = (weight[i]>0) ? kappahat[i]/weight[i] : 0;
    kappahat[i] += log_decay[i];
    distance.push_back(i+1);
  }
  return DecaySummary{distance,kappahat,weight};
}

DecayFit pointwise_log_decay_fit(const DecaySummary& dec) {
  //forbid increase past first diagonal
  std::vector<double> log_decay = dec.kappahat;
  unsigned N = log_decay.size();
  for (unsigned i=0; i<N; ++i) {
    if (i>=2 && log_decay[i] > log_decay[i-1]) log_decay[i] = log_decay[i-1];
  }
  //compute weighted mean
  double avg=0;
  double wsum=0;
  for (unsigned i=0; i<N; ++i) {
    avg += log_decay[i]*dec.weight[i];
    wsum += dec.weight[i];
  }
  avg = avg / wsum;
  //no smoothing, subtract average and return
  for (unsigned i=0; i<N; ++i) {
    log_decay[i] -= avg;
  }
  return DecayFit{log_decay,dec,-1};
}

DecayFit spline_log_decay_fit(const DecaySummary& dec, double tol_val, unsigned Kdiag, unsigned max_iter, double sigma) {
  //extract data
  const Eigen::Map<const Eigen::VectorXd> y(dec.kappahat.data(),dec.kappahat.size());
  const Eigen::Map<const Eigen::VectorXd> w(dec.weight.data(),dec.weight.size());
  const Eigen::VectorXd S = w.array().inverse().sqrt().matrix();
  //X: build spline base on log distance
  const Eigen::Map<const Eigen::ArrayXd> distance(dec.distance.data(), dec.distance.size());
  const Eigen::VectorXd log_distance = distance.log().matrix();
  const Eigen::SparseMatrix<double> X = generate_spline_base(log_distance, Kdiag);
  //D: build difference matrix
  const Eigen::SparseMatrix<double> D = second_order_difference_matrix(Kdiag);
  //C: build constraint matrix to forbid increase
  const Eigen::SparseMatrix<double> Cin = - first_order_difference_matrix(Kdiag);
  
  //iteratively fit GAM on decay and estimate lambda
  GeneralizedAdditiveModel gam(y,S,X,D,sigma);
  gam.set_inequality_constraints(Cin);
  gam.optimize(max_iter,tol_val);
  //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
  Eigen::VectorXd log_decay = gam.get_mean();
  const double lambda = gam.get_lambda();
  
  //compute weighted mean
  double avg = w.dot(log_decay)/w.sum();
  //subtract average and return
  log_decay.array() -= avg;
  return DecayFit{std::vector<double>(log_decay.data(),log_decay.data()+log_decay.rows()),dec,lambda};
}

//one IRLS iteration for log decay, with a poisson model
std::vector<double> step_log_decay(const FastSignalData& data, double tol_val) {
  //compute summary statistics for decay
  DecaySummary dec = get_decay_summary(data);
  //infer new decay from summaries
  //DecayFit fit = pointwise_log_decay_fit(dec);
  DecayFit fit = spline_log_decay_fit(dec, tol_val);
  /*Rcpp::Rcout << "distance log_decay kappahat weight\n";
  for (unsigned i=0; i<fit.log_decay.size(); ++i)
  Rcpp::Rcout << fit.dec.distance[i] << " " << fit.log_decay[i] << " "
              << fit.dec.kappahat[i] << " " << fit.dec.weight[i] << "\n";*/
  return fit.log_decay;
}

}
}

