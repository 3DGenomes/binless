#ifndef FAST_BIAS_HPP
#define FAST_BIAS_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "spline.hpp"
#include "gam.hpp"
#include "AnalyticalGAMLibrary.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

struct BiasGAMConfig {
  BiasGAMConfig(double tol_val, double constraint_every) : tol_val(tol_val), constraint_every(constraint_every) {}
  
  //parameters for bias calculation
  double bf_per_kb=5;
  unsigned max_iter=100;
  double sigma=1e-4;
  unsigned bins_per_bf=100;
  
  //default values will be overwritten
  double tol_val;      // tolerance on bias convergence
  unsigned constraint_every; // if > 0, the bias is forced to be centered every so many bases
  
};

class BiasGAMSettings {
  
public:
  template<typename FastData>
  BiasGAMSettings(const FastData& data, const BiasGAMConfig& conf) : conf_(conf) {
    //log_distance and bounds
    auto pos1_std = data.get_pos1();
    auto pos2_std = data.get_pos2();
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > pos1_data(pos1_std.data(),pos1_std.size());
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > pos2_data(pos2_std.data(),pos2_std.size());
    Eigen::VectorXd pos_data(pos1_data.rows()+pos2_data.rows());
    pos_data << pos1_data.cast<double>(), pos2_data.cast<double>();
    pos_min_ = pos_data.minCoeff();
    pos_max_ = pos_data.maxCoeff();
    //binner matrix
    K_ = conf_.bf_per_kb*(pos_max_-pos_min_)/1000.;
    const unsigned Nbins = K_ * conf_.bins_per_bf;
    const auto design = make_even_bins(pos_min_, pos_max_, Nbins);
    const bool drop = true; //drop unused bins
    auto binner1 = bin_data(pos1_data.cast<double>(), design, drop);
    auto binner2 = bin_data(pos2_data.cast<double>(), design, drop);
    binner_ = (  binner1 + binner2 );
    nbins_ = binner_.rows();
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    nobs_ = binner_ * nobs_data / 2.; // each obs is used twice in the binner matrix
    Rcpp::Rcout << "verif: sum(nobs_data)=" << nobs_data.sum() << " sum(nobs_)=" << nobs_.sum() << "\n";
    //compute mean position (currently, positions dont change within a bin but that might evolve)
    position_ =((  binner1*(pos1_data.cast<double>().array()*nobs_data.array()).matrix()
                 + binner2*(pos1_data.cast<double>().array()*nobs_data.array()).matrix() ).array() / (2*nobs_).array()).matrix();
    Rcpp::Rcout << "verif: min(position_)=" << position_.minCoeff() << " max(position_)=" << position_.maxCoeff();
    Rcpp::Rcout << " pos_min_=" << pos_min_ << " pos_max_=" << pos_max_ << "\n";
    Rcpp::Rcout << position_.head(5);
    Rcpp::Rcout << position_.tail(5);
    //X: design matrix
    X_ = generate_spline_base(position_, pos_min_, pos_max_, K_);
    //D: build difference matrix
    D_ = second_order_difference_matrix(K_);
    //C: build constraint matrix to center
    if (conf_.constraint_every > 0) {
      unsigned constraint_number = (pos_max_-pos_min_)/conf_.constraint_every;
      Rcpp::Rcout << "Using " << constraint_number << "=" << (pos_max_-pos_min_) << "/" << conf_.constraint_every << " constraints to model bias\n";
      const auto design2 = make_even_bins(pos_min_, pos_max_, constraint_number);
      Ceq_ = bin_data(position_, design2, drop);
    } else {
      Ceq_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(0,K_);
    }
  }
  
  double get_K() const { return K_; }
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
  const BiasGAMConfig& conf_;
  Eigen::VectorXd position_;
  double pos_min_, pos_max_;
  unsigned K_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd nobs_;
  Eigen::SparseMatrix<double> X_,D_,Ceq_; // design, difference and constraint matrices
};

struct BiasSummary {
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, phihat);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, weight);
};

struct BiasParams {
  BiasParams(const BiasGAMSettings& settings) : beta_(Eigen::VectorXd::Zero(settings.get_K())), lambda_(-1), mean_(0) {}

  BiasParams(const Rcpp::List& state) { set_state(state); }
  Rcpp::List get_state() const {
    return Rcpp::List::create(_["beta"]=get_beta(), _["lambda"]=get_lambda(), _["mean"]=get_mean());
  }
  void set_state(const Rcpp::List& state) {
    set_beta(Rcpp::as<Eigen::VectorXd>(state["beta"]));
    set_lambda(Rcpp::as<double>(state["lambda"]));
    set_mean(Rcpp::as<double>(state["mean"]));
  }
  
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, beta);
  BINLESS_GET_SET_DECL(double, double, lambda);
  BINLESS_GET_SET_DECL(double, double, mean);
};

class BiasGAMEstimator {
public:
  
  //here, initialize flat log bias
  template<typename FastData>
  BiasGAMEstimator(const FastData& data, const BiasGAMConfig& conf) :
   settings_(data,conf), summary_(), params_(settings_),
   gam_(settings_.get_X(), settings_.get_D(), settings_.get_sigma())
    { gam_.set_equality_constraints(settings_.get_Ceq()); }
  
  
  //compute group sums of a vector of the size of the input data into the bins formed for the bias calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return settings_.get_binner()*vec; }
  
  //get log bias along binned distances
  Eigen::VectorXd get_binned_estimate() const {
    return get_estimate() - Eigen::VectorXd::Constant(settings_.get_nbins(), params_.get_mean());
  }
  
  //get approximate log bias (bi + bj) along distances in original data (same approx as during fitting)
  Eigen::VectorXd get_data_estimate() const {
    return settings_.get_binner().transpose()*get_binned_estimate();
  }
  
  //provide a way to store and recall the state of the estimator
  Rcpp::List get_state() const { return params_.get_state(); }
  void set_state(const Rcpp::List& state) { params_.set_state(state); }
  
  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const std::vector<double>& log_expected, const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z);
  //perform spline fit of summary data
  void update_params();
  
private:
  //get X*beta
  Eigen::VectorXd get_estimate() const { return settings_.get_X() * params_.get_beta(); }
  
  Eigen::VectorXd get_beta() const { return params_.get_beta(); }
  void set_beta(const Eigen::VectorXd& beta) { params_.set_beta(beta); }
  
  double get_lambda() const { return params_.get_lambda(); }
  void set_lambda(double lambda) { params_.set_lambda(lambda); }
  
  //compute average log bias (weighted by nobs) in order to center it
  void center_estimate() {
    params_.set_mean(settings_.get_nobs().dot(settings_.get_X() * params_.get_beta())/settings_.get_nobs().sum());
  }
  
  const BiasGAMSettings settings_; // parameters for performing the binning, constant
  BiasSummary summary_; // transformed data, iteration-specific
  BiasParams params_; // resulting fit, iteration-specific
  GeneralizedAdditiveModel<AnalyticalGAMLibrary> gam_; //used to fit parameters
};

//typedef BiasGAMConfig BiasConfig;
//typedef BiasGAMSettings BiasSettings;
//typedef BiasGAMEstimator BiasEstimator;

}
}

#endif

