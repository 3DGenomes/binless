#ifndef FAST_DECAY_HPP
#define FAST_DECAY_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "macros.hpp" //BINLESS_*
#include "spline.hpp" //generate_spline_base
#include "QuadProgGAMLibrary.hpp"
#include "fast_estimator.hpp" //Estimator, Fitter, Summarizer and implementations
#include "Traits.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

template<>
struct Config<Decay,GAM> {
  Config(double tol_val, double free_decay) : tol_val(tol_val), free_decay(free_decay) {}
  
  //parameters for decay calculation
  unsigned K=50;
  unsigned max_iter=100;
  double sigma=1;
  unsigned bins_per_bf=100;
  
  //default values will be overwritten
  double tol_val;      // tolerance on decay convergence
  unsigned free_decay; // distance in bases until which decay is not forced to decrease
  
};

template<>
class SummarizerSettings<Decay,GAM> {
  
public:
  template<typename FastData>
  SummarizerSettings(const FastData& data, const Config<Decay,GAM>& conf) {
    //log_distance and bounds
    auto distance_std = data.get_distance();
    const Eigen::Map<const Eigen::VectorXd> distance_data(distance_std.data(),distance_std.size());
    Eigen::VectorXd log_distance_data = distance_data.array().log();
    log_dmin_ = log_distance_data.minCoeff();
    log_dmax_ = log_distance_data.maxCoeff();
    //binner matrix
    binner_ = bin_data_evenly(log_distance_data, conf.K*conf.bins_per_bf, true); //true to drop unused bins
    nbins_ = binner_.rows();
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    nobs_ = binner_ * nobs_data;
    //compute mean log distance
    log_distance_ = ((binner_ * (log_distance_data.array() * nobs_data.array()).matrix()).array() / nobs_.array()).matrix();
  }
  
  Eigen::VectorXd get_support() const { return log_distance_; }
  double get_support_min() const { return log_dmin_; }
  double get_support_max() const { return log_dmax_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  unsigned get_nbins() const { return nbins_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  
  BINLESS_FORBID_COPY(SummarizerSettings);
  
private:
  Eigen::VectorXd log_distance_;
  double log_dmin_, log_dmax_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd nobs_;
};

template<>
class FitterSettings<Decay,GAM> {
  
public:
  FitterSettings(const SummarizerSettings<Decay,GAM>& settings, const Config<Decay,GAM>& conf) :
      max_iter_(conf.max_iter), tol_val_(conf.tol_val), sigma_(conf.sigma), K_(conf.K), nbins_(settings.get_nbins()), nobs_(settings.get_nobs()) {
    auto log_distance = settings.get_support();
    auto log_dmin = settings.get_support_min();
    auto log_dmax = settings.get_support_max();
    //
    X_ = generate_spline_base(log_distance, log_dmin, log_dmax, get_K());
    //D: build difference matrix
    D_ = second_order_difference_matrix(get_K());
    //C: build constraint matrix to forbid increase
    unsigned free_first = get_K() * (std::log(conf.free_decay)-log_dmin)/(log_dmax-log_dmin);
    //Rcpp::Rcout << "Free decay in " << free_first << " out of " << conf_.K << " basis functions\n";
    Cin_ = decreasing_constraint(conf.K, free_first);
  }
  
  //this estimate is always centered after fitting
  bool is_centered() const { return true; }
  
  BINLESS_GET_CONSTREF_DECL(unsigned, max_iter);
  BINLESS_GET_CONSTREF_DECL(double, tol_val);
  BINLESS_GET_CONSTREF_DECL(double, sigma);
  BINLESS_GET_CONSTREF_DECL(double, K);
  BINLESS_GET_CONSTREF_DECL(unsigned, nbins);
  BINLESS_GET_CONSTREF_DECL(Eigen::VectorXd, nobs);
  
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, X);
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, D);
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, Cin);
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, Ceq); //not used
  
  BINLESS_FORBID_COPY(FitterSettings);
};

template<>
struct FitterTraits<Decay,GAM> {
  typedef QuadProgGAMLibrary library;
  static const bool has_inequality_constraints = true;
  static const bool has_equality_constraints = false;
};


typedef Config<Decay,GAM> DecayConfig;
typedef Estimator<Decay,GAM> DecayEstimator;

}
}

#endif

