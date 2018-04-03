#ifndef FAST_BIAS_MEAN_HPP
#define FAST_BIAS_MEAN_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "spline.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

struct BiasMeanConfig {
  BiasMeanConfig(unsigned nbins) : nbins(nbins) {}
  unsigned nbins;
};

class BiasMeanSettings {
  
public:
  template<typename FastData>
  BiasMeanSettings(const FastData& data, const BiasMeanConfig& conf) : conf_(conf) {
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
    const auto design = make_even_bins(pos_min_, pos_max_, conf_.nbins);
    const bool drop = true; //drop unused bins
    auto binner1 = bin_data(pos1_data.cast<double>(), design, drop);
    auto binner2 = bin_data(pos2_data.cast<double>(), design, drop);
    binner_ = (  binner1 + binner2 );
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
  }
  
  unsigned get_nbins() const { return conf_.nbins; }
  
  Eigen::VectorXd get_position() const { return position_; }
  double get_pos_min() const { return pos_min_; }
  double get_pos_max() const { return pos_max_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  
private:
  const BiasMeanConfig& conf_;
  Eigen::VectorXd position_;
  double pos_min_, pos_max_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  Eigen::VectorXd nobs_;
};

struct BiasMeanSummary {
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, etahat);
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, weight);
};

struct BiasMeanParams {
  BiasMeanParams(const BiasMeanSettings& settings) : estimate_(Eigen::VectorXd::Zero(settings.get_nbins())), mean_(0) {}
  
  BINLESS_GET_SET_DECL(Eigen::VectorXd, const Eigen::VectorXd&, estimate);
  BINLESS_GET_SET_DECL(double, double, mean);
};

class BiasMeanEstimator {
public:
  
  //here, initialize flat log bias
  template<typename FastData>
  BiasMeanEstimator(const FastData& data, const BiasMeanConfig& conf) :
   settings_(data,conf), summary_(), params_(settings_) {}
  
  
  //compute group sums of a vector of the size of the input data into the bins formed for the bias calculation
  Eigen::VectorXd summarize(const Eigen::VectorXd& vec) const { return settings_.get_binner()*vec; }
  
  //set log bias
  void set_estimate(const Eigen::VectorXd& log_bias) { params_.set_estimate(log_bias); }
  
  //get log bias along binned distances
  Eigen::VectorXd get_binned_estimate() const {
    return params_.get_estimate() - Eigen::VectorXd::Constant(settings_.get_nbins(), params_.get_mean());
  }
  
  //get approximate log bias (bi + bj) along distances in original data (same approx as during fitting)
  Eigen::VectorXd get_data_estimate() const {
    return settings_.get_binner().transpose()*get_binned_estimate();
  }
  
  //initial guess of IRLS weights using poisson model
  void set_poisson_lsq_summary(const std::vector<double>& log_expected, const FastSignalData& data, double pseudocount=0.01);
  //incremental update of IRLS weights
  void update_summary(const ResidualsPair& z);
  //perform spline fit of summary data
  void update_params();
  
private:
  //compute average log bias (weighted by nobs) in order to center it
  void center_estimate() {
    params_.set_mean(settings_.get_nobs().dot(params_.get_estimate())/settings_.get_nobs().sum());
  }
  
  const BiasMeanSettings settings_; // parameters for performing the binning, constant
  BiasMeanSummary summary_; // transformed data, iteration-specific
  BiasMeanParams params_; // resulting fit, iteration-specific
};

typedef BiasMeanConfig BiasConfig;
typedef BiasMeanSettings BiasSettings;
typedef BiasMeanEstimator BiasEstimator;

}
}

#endif

