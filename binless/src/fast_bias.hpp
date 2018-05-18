#ifndef FAST_BIAS_HPP
#define FAST_BIAS_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "FastData.hpp"
#include "util.hpp" //bin_data_evenly
#include "macros.hpp"
#include "spline.hpp"
#include "gam.hpp"
#include "AnalyticalGAMLibrary.hpp"
#include "fast_estimator.hpp"
#include "Traits.hpp"

namespace binless {
namespace fast {

struct ResidualsPair;

template<>
struct Config<Bias,GAM> {
  Config(double tol_val, double constraint_every) : tol_val(tol_val), constraint_every(constraint_every) {}
  
  //parameters for bias calculation
  double bf_per_kb=5;
  unsigned max_iter=100;
  double sigma=1e-4;
  unsigned bins_per_bf=100;
  
  //default values will be overwritten
  double tol_val;      // tolerance on bias convergence
  unsigned constraint_every; // if > 0, the bias is forced to be centered every so many bases
  
};

//enable debug
template<>
struct SummarizerTraits<Bias,GAM> {
  //print debug info
  static const bool debug = false;
};

template<>
class SummarizerSettingsImpl<Bias,GAM> : public SummarizerSettings {
  
public:
  template<typename FastData>
  SummarizerSettingsImpl(const FastData& data, const Config<Bias,GAM>& conf) {
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
    const unsigned K = conf.bf_per_kb*(get_support_max()-get_support_min())/1000.;
    const unsigned Nbins = K * conf.bins_per_bf;
    const auto design = make_even_bins(get_support_min(), get_support_max(), Nbins);
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
    set_support(   ((  binner1*(pos1_data.cast<double>().array()*nobs_data.array()).matrix()
                     + binner2*(pos1_data.cast<double>().array()*nobs_data.array()).matrix() ).array() / (2*get_nobs()).array()).matrix() );
  }
};

template<>
class FitterSettingsImpl<Bias,GAM> : public FitterSettings<GAM> {
  
public:
  FitterSettingsImpl(const SummarizerSettings& settings, const Config<Bias,GAM>& conf) :
      FitterSettings<GAM>(conf.max_iter, conf.tol_val, conf.sigma,
                          conf.bf_per_kb*(settings.get_support_max()-settings.get_support_min())/1000., // K
                          settings.get_nbins(), settings.get_nobs()) {
    auto position = settings.get_support();
    auto pos_min = settings.get_support_min();
    auto pos_max = settings.get_support_max();
    //X: design matrix
    set_X( generate_spline_base(position, pos_min, pos_max, get_K()) );
    //D: build difference matrix
    set_D( second_order_difference_matrix(get_K()) );
    //C: build constraint matrix to center
    if (conf.constraint_every > 0) {
      unsigned constraint_number = (pos_max-pos_min)/conf.constraint_every;
      Rcpp::Rcout << "Using " << constraint_number << "=" << (pos_max-pos_min) << "/" << conf.constraint_every << " constraints to model bias\n";
      const auto design2 = make_even_bins(pos_min, pos_max, constraint_number);
      const bool drop = true;
      set_Ceq( bin_data(position, design2, drop) );
    } else {
      set_Ceq( Eigen::SparseMatrix<double, Eigen::RowMajor>(0,get_K()) );
    }
  }
};

template<>
struct FitterTraits<Bias,GAM> {
  typedef AnalyticalGAMLibrary library;
  static const bool has_inequality_constraints = false;
  static const bool has_equality_constraints = true;
  //print debug info
  static const bool debug = false;
};

//typedef Config<Bias,GAM> BiasConfig;
//typedef Estimator<Bias,GAM> BiasEstimator;

}
}

#endif

