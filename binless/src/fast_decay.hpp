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

//enable debug
template<>
struct SummarizerTraits<Decay,GAM> {
  //print debug info
  static const bool debug = false;
};

template<>
class SummarizerSettingsImpl<Decay,GAM> : public SummarizerSettings {
  
public:
  template<typename FastData>
  SummarizerSettingsImpl(const FastData& data, const Config<Decay,GAM>& conf) {
    //log_distance and bounds
    auto distance_std = data.get_distance();
    const Eigen::Map<const Eigen::VectorXd> distance_data(distance_std.data(),distance_std.size());
    Eigen::VectorXd log_distance_data = distance_data.array().log();
    set_support_min(log_distance_data.minCoeff());
    set_support_max(log_distance_data.maxCoeff());
    //binner matrix
    set_binner(bin_data_evenly(log_distance_data, conf.K*conf.bins_per_bf, true)); //true to drop unused bins
    set_nbins(get_binner().rows());
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    set_nobs(get_binner() * nobs_data);
    //compute mean log distance
    //in theory, should be weighted by nobs, but that causes problems when nobs=0
    //counter diagonals with nobs=0 cannot yet be removed because signal is not ported to new class hierarchy
    set_support( (get_binner() * log_distance_data).cwiseQuotient(
                        get_binner() * Eigen::VectorXd::Ones(log_distance_data.rows())) );
  }
};

template<>
class FitterSettingsImpl<Decay,GAM> : public FitterSettings<GAM> {
  
public:
  FitterSettingsImpl(const SummarizerSettings& settings, const Config<Decay,GAM>& conf) :
      FitterSettings<GAM>(conf.max_iter, conf.tol_val, conf.sigma, conf.K, settings.get_nbins(), settings.get_nobs()) {
    auto log_distance = settings.get_support();
    auto log_dmin = settings.get_support_min();
    auto log_dmax = settings.get_support_max();
    //
    set_X( generate_spline_base(log_distance, log_dmin, log_dmax, get_K()) );
    //D: build difference matrix
    set_D( second_order_difference_matrix(get_K()) );
    //C: build constraint matrix to forbid increase
    unsigned free_first = get_K() * (std::log(conf.free_decay)-log_dmin)/(log_dmax-log_dmin);
    //Rcpp::Rcout << "Free decay in " << free_first << " out of " << conf_.K << " basis functions\n";
    set_Cin( decreasing_constraint(get_K(), free_first) );
  }
};

template<>
struct FitterTraits<Decay,GAM> {
  typedef QuadProgGAMLibrary library;
  static const bool has_inequality_constraints = true;
  static const bool has_equality_constraints = false;
  //print debug info
  static const bool debug = false;
};


typedef Config<Decay,GAM> DecayConfig;
typedef Estimator<Decay,GAM> DecayEstimator;

}
}

#endif

