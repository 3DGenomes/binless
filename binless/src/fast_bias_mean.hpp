#ifndef FAST_BIAS_MEAN_HPP
#define FAST_BIAS_MEAN_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "macros.hpp"
#include "spline.hpp"
#include "fast_estimator.hpp"

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
  
  //this estimate is never centered after fitting (internal centering)
  bool is_centered() const { return false; }
  
  Eigen::VectorXd get_position() const { return position_; }
  double get_pos_min() const { return pos_min_; }
  double get_pos_max() const { return pos_max_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  
  BINLESS_FORBID_COPY(BiasMeanSettings);
  
private:
  const BiasMeanConfig& conf_;
  Eigen::VectorXd position_;
  double pos_min_, pos_max_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  Eigen::VectorXd nobs_;
};

class BiasMeanImpl {
public:
  template<typename FastData>
  BiasMeanImpl(const FastData& data, const BiasMeanConfig& conf) : 
    settings_(data, conf), summary_(), params_(settings_) {}
  
  void update_params();
  
  //get X*beta
  Eigen::VectorXd get_estimate() const { return params_.get_estimate(); }
  
  BINLESS_GET_CONSTREF_DECL(BiasMeanSettings, settings);
  BINLESS_GET_REF_DECL(Summary, summary);
  BINLESS_GET_REF_DECL(MeanParams, params);
  
private:
  void set_estimate(const Eigen::VectorXd& estimate) { params_.set_estimate(estimate); }
  
};

typedef BiasMeanConfig BiasConfig;
typedef BiasMeanSettings BiasSettings;
typedef Estimator<BiasMeanImpl> BiasEstimator;

}
}

#endif

