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

template<>
class SummarizerSettings<Bias,Mean> {
  
public:
  template<typename FastData>
  SummarizerSettings(const FastData& data, const Config<Bias,Mean>& conf) {
    //log_distance and bounds
    auto pos1_std = data.get_pos1();
    auto pos2_std = data.get_pos2();
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > pos1_data(pos1_std.data(),pos1_std.size());
    const Eigen::Map<const Eigen::Matrix<unsigned,Eigen::Dynamic,1> > pos2_data(pos2_std.data(),pos2_std.size());
    Eigen::VectorXd pos_data(pos1_data.rows()+pos2_data.rows());
    pos_data << pos1_data.cast<double>(), pos2_data.cast<double>();
    support_min_ = pos_data.minCoeff();
    support_max_ = pos_data.maxCoeff();
    //binner matrix
    const auto design = make_even_bins(get_support_min(), get_support_max(), conf.nbins);
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
    //compute mean position (currently, positions dont change within a bin but that might evolve)
    support_ =((  binner1*(pos1_data.cast<double>().array()*nobs_data.array()).matrix()
                     + binner2*(pos1_data.cast<double>().array()*nobs_data.array()).matrix() ).array() / (2*nobs_).array()).matrix();
  }
  
  Eigen::VectorXd get_support() const { return support_; }
  double get_support_min() const { return support_min_; }
  double get_support_max() const { return support_max_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  unsigned get_nbins() const { return nbins_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  
  BINLESS_FORBID_COPY(SummarizerSettings);
  
private:
  Eigen::VectorXd support_;
  double support_min_, support_max_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd nobs_;
};

template<>
class FitterSettings<Bias,Mean> {
  
public:
  FitterSettings(const SummarizerSettings<Bias,Mean>& settings, const Config<Bias,Mean>&) :
    nbins_(settings.get_nbins()), nobs_(settings.get_nobs()) {}
    
  //this estimate is never centered after fitting (internal centering)
  bool is_centered() const { return false; }
  
  BINLESS_GET_CONSTREF_DECL(unsigned, nbins);
  BINLESS_GET_CONSTREF_DECL(Eigen::VectorXd, nobs);
  BINLESS_FORBID_COPY(FitterSettings);
};

template<>
struct FitterTraits<Bias,Mean> {};

typedef Config<Bias,Mean> BiasConfig;
typedef Estimator<Bias,Mean> BiasEstimator;

}
}

#endif

