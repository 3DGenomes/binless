#ifndef FAST_EXPOSURE_HPP
#define FAST_EXPOSURE_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "macros.hpp"
#include "spline.hpp"
#include "fast_estimator.hpp"
#include "Traits.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

template<>
struct Config<Exposure,Mean> {};


//enable debug
template<>
struct SummarizerTraits<Exposure,Mean> {
  //print debug info
  static const bool debug = false;
};

template<>
class SummarizerSettingsImpl<Exposure,Mean> : public SummarizerSettings {
  
public:
  template<typename FastData>
  SummarizerSettingsImpl(const FastData& data, const Config<Exposure,Mean>&) {
    //ids for each dataset
    auto name_std = data.get_name();
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > name_data(name_std.data(),name_std.size());
    set_support_min(name_data.minCoeff());
    set_support_max(name_data.maxCoeff());
    //binner matrix
    set_binner(bin_data_evenly(name_data.cast<double>(), get_support_max()-get_support_min()+1, true)); //true to drop unused bins
    set_nbins(get_binner().rows());
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    set_nobs(get_binner() * nobs_data);
    //compute mean log distance
    set_support( ((get_binner() * (name_data.cast<double>().array() * nobs_data.array()).matrix()).array() / get_nobs().array()).matrix() );
    Rcpp::Rcout << "support: " << get_support().transpose() << "\n";
  }
};

template<>
class FitterSettingsImpl<Exposure,Mean> : public FitterSettings<Mean> {
  
public:
  FitterSettingsImpl(const SummarizerSettings& settings, const Config<Exposure,Mean>&) :
    FitterSettings<Mean>(settings.get_nbins(), settings.get_nobs()) {}
};

template<>
struct FitterTraits<Exposure,Mean> {
  //do not center the estimate!
  static const bool center = false;
  //do not cap the data, use exact means
  static const bool cap = false;
  //print debug info
  static const bool debug = false;
};

typedef Config<Exposure,Mean> ExposureConfig;
typedef Estimator<Exposure,Mean> ExposureEstimator;

}
}

#endif

