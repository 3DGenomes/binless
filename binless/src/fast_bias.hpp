#ifndef FAST_BIAS_HPP
#define FAST_BIAS_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "GFLLibrary.hpp"
#include "util.hpp" //bin_data_evenly
#include "spline.hpp"
#include "gam.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

struct BiasConfig {
  BiasConfig(double tol_val, double constraint_every) : tol_val(tol_val), constraint_every(constraint_every) {}
  
  //parameters for bias calculation
  double bf_per_kb=5;
  unsigned max_iter=100;
  double sigma=1e-4;
  unsigned bins_per_bf=100;
  
  //default values will be overwritten
  double tol_val;      // tolerance on bias convergence
  unsigned constraint_every; // if constrained, the bias is forced to be centered every so many bases
  
};

class BiasSettings {
  
public:
  template<typename FastData>
  BiasSettings(const FastData& data, const BiasConfig& conf) : conf_(conf) {
    //log_distance and bounds
    auto pos1_std = data.get_pos1();
    auto pos2_std = data.get_pos2();
    const Eigen::Map<const Eigen::VectorXd> pos1_data(pos1_std.data(),pos1_std.size());
    const Eigen::Map<const Eigen::VectorXd> pos2_data(pos2_std.data(),pos2_std.size());
    Eigen::VectorXd pos_data(pos1_data.rows()+pos2_data.rows());
    pos_data << pos1_data, pos2_data;
    pos_min_ = pos_data.minCoeff();
    pos_max_ = pos_data.maxCoeff();
    //binner matrix
    Krow_ = conf_.bf_per_kb*(pos_max_-pos_min_)/1000.;
    const unsigned Nbins = Krow_ * conf_.bins_per_bf;
    const auto design = make_even_bins(pos_min_, pos_max_, Nbins);
    const bool drop = true; //drop unused bins
    binner_ = (  bin_data(pos1_data, design, drop)
               + bin_data(pos2_data, design, drop) )/2.; //divide by 2 because we use the data twice
    nbins_ = binner_.rows();
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    nobs_ = binner_ * nobs_data;
    //compute mean position (currently, positions dont change within a bin but that might evolve)
    position_ = ((binner_ * (pos_data.array() * nobs_data.array()).matrix()).array() / nobs_.array()).matrix();
    //X: design matrix
    X_ = generate_spline_base(position_, pos_min_, pos_max_, Krow_);
    //D: build difference matrix
    D_ = second_order_difference_matrix(Krow_);
    //C: build constraint matrix to center
    unsigned constraint_number = (pos_max_-pos_min_)/conf_.constraint_every;
    Rcpp::Rcout << "Using " << constraint_number << "=" << (pos_max_-pos_min_) << "/" << conf_.constraint_every << " constraints to model bias\n";
    const auto design2 = make_even_bins(pos_min_, pos_max_, constraint_number);
    Ceq_ = bin_data(position_, design2, drop);
  }
  
  double get_Krow() const { return Krow_; }
  unsigned get_nbins() const { return nbins_; }
  unsigned get_max_iter() const { return conf_.max_iter; }
  double get_tol_val() const { return conf_.tol_val; }
  double get_sigma() const { return conf_.sigma; }
  
  Eigen::VectorXd get_position() const { return position_; }
  double get_pos_min() const { return pos_min_; }
  double get_pos_max() const { return pos_max_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  Eigen::SparseMatrix<double> get_X() const { return X_; }
  Eigen::SparseMatrix<double> get_D() const { return D_; }
  Eigen::SparseMatrix<double> get_Ceq() const { return Ceq_; }
  
private:
  const BiasConfig& conf_;
  Eigen::VectorXd position_;
  double pos_min_, pos_max_;
  unsigned Krow_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd nobs_;
  Eigen::SparseMatrix<double> X_,D_,Ceq_; // design, difference and constraint matrices
};

struct BiasSummary {
  Eigen::VectorXd etahat, weight;
};

struct BiasParams {
  BiasParams(const BiasSettings& settings) : beta(Eigen::VectorXd::Zero(settings.get_Krow())), lambda(-1), mean(0) {}
  
  Eigen::VectorXd beta;
  double lambda, mean;
};

class BiasEstimator {
public:
  
  //here, initialize flat log bias
  template<typename FastData>
  BiasEstimator(const FastData& data, const BiasConfig& conf) :
   settings_(data,conf), summary_(), params_(settings_),
   gam_(settings_.get_X(), settings_.get_D(), settings_.get_sigma())
    { gam_.set_equality_constraints(settings_.get_Ceq()); }
  
  
  //compute group sums of a vector of the size of the input data into the bins formed for the bias calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return settings_.get_binner()*vec; }
  
  Eigen::VectorXd get_beta() const { return params_.beta; }
  void set_beta(const Eigen::VectorXd& beta) { params_.beta = beta; }
  
  //get log bias along binned distances
  Eigen::VectorXd get_binned_estimate() const {
    return settings_.get_X() * params_.beta - Eigen::VectorXd::Constant(settings_.get_nbins(), params_.mean);
  }
  
  //get approximate log bias along distances in original data (same approx as during fitting)
  Eigen::VectorXd get_data_estimate() const {
    return settings_.get_binner().transpose()*get_binned_estimate();
  }
  
  double get_lambda() const { return params_.lambda; }
  void set_lambda(double lambda) { params_.lambda = lambda; }

  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z);
  //perform spline fit of summary data
  void update_params();
  //one complete IRLS iteration for log bias
  void step_irls(const ResidualsPair& z) {
    update_summary(z);
    update_params();
  }
  
  
private:
  //compute average log bias (weighted by nobs) in order to center it
  void center_estimate() {
    params_.mean = settings_.get_nobs().dot(settings_.get_X() * params_.beta)/settings_.get_nobs().sum();
  }
  
  const BiasSettings settings_; // parameters for performing the binning, constant
  BiasSummary summary_; // transformed data, iteration-specific
  BiasParams params_; // resulting fit, iteration-specific
  GeneralizedAdditiveModel gam_; //used to fit parameters
};

}
}

#endif

