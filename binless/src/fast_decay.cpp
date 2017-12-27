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
#include "fast_expected.hpp"

namespace binless {
namespace fast {


//here, log decay is log ( sum_i observed / sum_i expected ) with i summed over counter diagonals
DecaySummary compute_poisson_lsq_log_decay(const FastSignalData& data, const DecayEstimate& dec, const DecaySchedule& schedule) {
  auto distance_std = data.get_distance();
  const Eigen::Map<const Eigen::VectorXd> distance(distance_std.data(),distance_std.size());
  Eigen::ArrayXd log_distance = distance.array().log();
  const auto bin_mat = bin_data_evenly(log_distance, schedule.bins_per_bf*schedule.Kdiag, true); //true to drop unused bins
  //compute ncounts
  Eigen::VectorXd ncounts = bin_mat * Eigen::VectorXd::Ones(log_distance.size());
  //compute mean distance
  Eigen::VectorXd mean_distance = ((bin_mat * log_distance.matrix()).array() / ncounts.array()).exp().matrix();
  //compute observed data
  auto observed_std = data.get_observed();
  Eigen::VectorXd observed = Eigen::VectorXd::Zero(observed_std.size());
  for (unsigned i=0; i<observed.rows(); ++i) observed(i) = observed_std[i]; // cast to double
  Eigen::VectorXd sum_obs = bin_mat * observed;
  //compute expected data
  auto log_expected_std = get_log_expected(data, dec);
  const Eigen::Map<const Eigen::VectorXd> log_expected(log_expected_std.data(),log_expected_std.size());
  Eigen::VectorXd sum_exp = bin_mat * (log_expected.array().exp().matrix());
  //compute log_decay
  for (unsigned i=0; i<data.get_nbins(); ++i) {
    if (sum_obs(i)==0) {
      Rcpp::Rcout << "counter diag " << i << " is zero!\n";
      Rcpp::stop(" Aborting...");
    }
  }
  Eigen::VectorXd log_decay = (sum_obs.array()/sum_exp.array()).log().matrix();
  //center log_decay
  double avg = log_decay.sum()/log_decay.rows();
  log_decay = (log_decay.array() - avg).matrix();
  //compute weight
  Eigen::VectorXd weight = log_decay.array().exp().matrix();
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance kappahat weight ncounts\n";
  Rcpp::Rcout << (Eigen::MatrixXd(mean_distance.rows(),4) << mean_distance, log_decay,
                  weight, ncounts).finished();*/
  return DecaySummary{mean_distance,log_decay,weight,ncounts};
}

DecaySummary get_decay_summary(const FastSignalData& data, const DecayEstimate& dec, const DecaySchedule& schedule) {
  //get residuals
  ResidualsPair z = get_poisson_residuals(data, dec);
  //build bins to agglomerate them
  auto distance_std = data.get_distance();
  const Eigen::Map<const Eigen::VectorXd> distance(distance_std.data(),distance_std.size());
  Eigen::ArrayXd log_distance = distance.array().log();
  const auto bin_mat = bin_data_evenly(log_distance, schedule.bins_per_bf*schedule.Kdiag, true); //true to drop unused bins
  //compute ncounts
  Eigen::VectorXd ncounts = bin_mat * Eigen::VectorXd::Ones(log_distance.size());
  //compute mean distance
  Eigen::VectorXd mean_distance = ((bin_mat * log_distance.matrix()).array() / ncounts.array()).exp().matrix();
  //compute weight
  const Eigen::Map<const Eigen::VectorXd> weights(z.weights.data(),z.weights.size());
  Eigen::VectorXd mean_weight = bin_mat * weights;
  //compute kappahat
  const Eigen::Map<const Eigen::VectorXd> residuals(z.residuals.data(),z.residuals.size());
  Eigen::VectorXd kappahat = bin_mat * (residuals.array() * weights.array()).matrix();
  for (unsigned i=0; i<kappahat.rows(); ++i) {
    if (mean_weight(i)>0) {
      kappahat(i) =  kappahat(i)/mean_weight(i);
    } else {
      kappahat(i) = 0;
    }
  }
  Eigen::VectorXd log_decay = dec.get_log_decay(mean_distance.array().log());
  kappahat += log_decay;
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance kappahat weight ncounts log_decay\n";
  Rcpp::Rcout << mean_distance.rows() << " " << kappahat.rows() << " " << mean_weight.rows() << " " << ncounts.rows() << "\n";
  Rcpp::Rcout << (Eigen::MatrixXd(mean_distance.rows(),5) << mean_distance, kappahat,
                  mean_weight, ncounts, log_decay).finished();*/
  return DecaySummary{mean_distance,kappahat,mean_weight,ncounts};
}

void spline_log_decay_fit(const DecaySummary& summary, DecayEstimate& dec, double tol_val, const DecaySchedule& schedule) {
  //extract data
  const Eigen::VectorXd y(summary.kappahat);
  const Eigen::VectorXd w(summary.weight);
  const Eigen::VectorXd S = w.array().inverse().sqrt().matrix();
  //X: build spline base on log distance
  const Eigen::ArrayXd distance(summary.distance.array());
  const Eigen::VectorXd log_distance = distance.log().matrix();
  const Eigen::SparseMatrix<double> X = dec.get_design(log_distance);
  //D: build difference matrix
  const Eigen::SparseMatrix<double> D = second_order_difference_matrix(schedule.Kdiag);
  //C: build constraint matrix to forbid increase
  const Eigen::SparseMatrix<double> Cin = - first_order_difference_matrix(schedule.Kdiag);
  
  //iteratively fit GAM on decay and estimate lambda
  GeneralizedAdditiveModel gam(y,S,X,D,schedule.sigma);
  gam.set_inequality_constraints(Cin);
  gam.optimize(schedule.max_iter,tol_val);
  //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
  dec.set_beta_diag(gam.get_beta());
  dec.center_log_decay(log_distance, summary.ncounts);
  dec.set_lambda_diag(gam.get_lambda());
  dec.set_summary(summary);
  /*Rcpp::Rcout << "spline_log_decay_fit\n";
  Rcpp::Rcout << "distance kappahat weight ncounts log_decay\n";
  Rcpp::Rcout << (Eigen::MatrixXd(summary.distance.rows(),5) << summary.distance, summary.kappahat,
                  summary.weight, summary.ncounts, dec.get_log_decay(summary.distance.array().log().matrix())).finished();*/
}

//one IRLS iteration for log decay, with a poisson model
void step_log_decay(const FastSignalData& data, DecayEstimate& dec, double tol_val) {
  //compute summary statistics for decay
  DecaySchedule schedule;
  DecaySummary summary = get_decay_summary(data, dec, schedule);
  //infer new decay from summaries
  spline_log_decay_fit(summary, dec, tol_val, schedule);
}

}
}

