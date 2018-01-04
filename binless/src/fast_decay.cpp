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
void Decay::set_poisson_lsq_summary(const FastSignalData& data) {
  //compute ncounts
  Eigen::VectorXd ncounts = schedule_.get_ncounts();
  //compute mean distance
  Eigen::VectorXd mean_distance = schedule_.get_log_distance().array().exp().matrix();
  //compute observed data
  auto observed_std = data.get_observed();
  Eigen::VectorXd observed = Eigen::VectorXd::Zero(observed_std.size());
  for (unsigned i=0; i<observed.rows(); ++i) observed(i) = observed_std[i]; // cast to double
  Eigen::VectorXd sum_obs = summarize(observed);
  //compute expected data
  auto log_expected_std = get_log_expected(data, *this);
  const Eigen::Map<const Eigen::VectorXd> log_expected(log_expected_std.data(),log_expected_std.size());
  Eigen::VectorXd sum_exp = summarize(log_expected.array().exp().matrix());
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
  summary_.kappahat = log_decay;
  summary_.weight = weight;
}

void Decay::update_summary(const FastSignalData& data) {
  //get residuals
  ResidualsPair z = get_poisson_residuals(data, *this);
  //compute weight
  const Eigen::Map<const Eigen::VectorXd> weights(z.weights.data(),z.weights.size());
  Eigen::VectorXd weight_sum = summarize(weights);
  //compute kappahat
  const Eigen::Map<const Eigen::VectorXd> residuals(z.residuals.data(),z.residuals.size());
  Eigen::VectorXd kappahat = summarize(residuals.array() * weights.array()).matrix();
  for (unsigned i=0; i<kappahat.rows(); ++i) {
    if (weight_sum(i)>0) {
      kappahat(i) =  kappahat(i)/weight_sum(i);
    } else {
      Rcpp::Rcout << "weight sum is zero for distance bin " << i << "\n";
      Rcpp::stop("aborting");
      kappahat(i) = 0;
    }
  }
  Eigen::VectorXd log_decay = get_binned_log_decay();
  kappahat += log_decay;
  /*Rcpp::Rcout << "BEFORE\n";
  Rcpp::Rcout << "distance kappahat weight ncounts log_decay\n";
  Rcpp::Rcout << (Eigen::MatrixXd(kappahat.rows(),5) << schedule_.get_log_distance().array().exp().matrix(), kappahat,
                  weight_sum, log_decay).finished();*/
  summary_.kappahat = kappahat;
  summary_.weight = weight_sum;
}

void Decay::update_params() {
  //extract data
  const Eigen::VectorXd y(summary_.kappahat);
  const Eigen::VectorXd w(summary_.weight);
  const Eigen::VectorXd S = w.array().inverse().sqrt().matrix();
  const Eigen::SparseMatrix<double> X = schedule_.get_X();
  const Eigen::SparseMatrix<double> D = schedule_.get_D();
  const Eigen::SparseMatrix<double> Cin = schedule_.get_Cin();
  
  //iteratively fit GAM on decay and estimate lambda
  GeneralizedAdditiveModel gam(y,S,X,D,schedule_.get_sigma());
  gam.set_inequality_constraints(Cin);
  gam.optimize(schedule_.get_max_iter(),schedule_.get_tol_val());
  //Rcpp::Rcout << "gam converged: " << gam.has_converged() << "\n";
  set_lambda_diag(gam.get_lambda());
  set_beta_diag(gam.get_beta());
  //center log decay
  double avg = schedule_.get_ncounts().dot(get_binned_log_decay())/schedule_.get_ncounts().sum();
  params_.mean += avg;
  /*Rcpp::Rcout << "spline_log_decay_fit\n";
  Rcpp::Rcout << "distance kappahat weight ncounts log_decay\n";
  Rcpp::Rcout << (Eigen::MatrixXd(schedule_.get_nbins(),5) << schedule.get_log_distance().array().exp().matrix(),
                  summary_.kappahat, summary_.weight, schedule_.get_ncounts(), get_binned_log_decay()).finished();*/
}

}
}

