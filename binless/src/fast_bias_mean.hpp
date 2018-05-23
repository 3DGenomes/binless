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
#include "Traits.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

template<>
struct Config<Bias,Mean> {
  Config(unsigned nbins) : nbins(nbins) {}
  unsigned nbins;
};

//enable debug
template<>
struct SummarizerTraits<Bias,Mean> {
  //print debug info
  static const bool debug = false;
};

template<>
class SummarizerSettingsImpl<Bias,Mean> : public SummarizerSettings {
  
public:
  template<typename FastData>
  SummarizerSettingsImpl(const FastData& data, const Config<Bias,Mean>& conf) {
    //log_distance and bounds
    auto pos1_std = data.get_pos1();
    auto pos2_std = data.get_pos2();
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > pos1_data(pos1_std.data(),pos1_std.size());
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > pos2_data(pos2_std.data(),pos2_std.size());
    Eigen::VectorXd pos_data(pos1_data.rows()+pos2_data.rows());
    pos_data << pos1_data.cast<double>(), pos2_data.cast<double>();
    set_support_min( pos_data.minCoeff() );
    set_support_max( pos_data.maxCoeff() );
    //binner matrix
    const auto design = make_even_bins(get_support_min(), get_support_max(), conf.nbins);
    const bool drop = true; //drop unused bins
    auto binner1 = bin_data(pos1_data.cast<double>(), design, drop);
    auto binner2 = bin_data(pos2_data.cast<double>(), design, drop);
    set_binner( (  binner1 + binner2 ) );
    set_nbins( get_binner().rows() );
    //nobs
    auto nobs_std = data.get_nobs();
    Eigen::VectorXd nobs_data = Eigen::VectorXd::Zero(nobs_std.size());
    for (unsigned i=0; i<nobs_data.rows(); ++i) nobs_data(i) = nobs_std[i]; // cast to double
    set_nobs( get_binner() * nobs_data / 2. ); // each obs is used twice in the binner matrix
    //compute mean position (currently, positions dont change within a bin but that might evolve)
    //in theory, should be weighted by nobs, but that causes problems when nobs=0
    //columns with nobs=0 cannot yet be removed because signal is not ported to new class hierarchy
    Eigen::VectorXd support = (   (binner1*pos1_data.cast<double>()).cwiseQuotient(
                                                    binner1*Eigen::VectorXd::Ones(pos1_data.rows()))
                                + (binner2*pos2_data.cast<double>()).cwiseQuotient(
                                                    binner2*Eigen::VectorXd::Ones(pos2_data.rows()))
                               )/2.;
    set_support(support);
  }
};

template<>
class FitterSettingsImpl<Bias,Mean> : public FitterSettings<Mean> {
  
public:
  FitterSettingsImpl(const SummarizerSettings& settings, const Config<Bias,Mean>&) :
    FitterSettings<Mean>(settings.get_nbins(), settings.get_nobs()) {}
};

template<>
struct FitterTraits<Bias,Mean> {
  //center the estimate
  static const bool center = true;
  //cap the data at 3SD
  static const bool cap = true;
  //print debug info
  static const bool debug = false;
};

typedef Config<Bias,Mean> BiasConfig;
typedef Estimator<Bias,Mean> BiasEstimator;

}
}

#endif

