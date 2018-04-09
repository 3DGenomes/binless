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

template<>
class SummarizerSettings<Bias,GAM> {
  
public:
  template<typename FastData>
  SummarizerSettings(const FastData& data, const Config<Bias,GAM>& conf) {
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
    const unsigned K = conf.bf_per_kb*(pos_max_-pos_min_)/1000.;
    const unsigned Nbins = K * conf.bins_per_bf;
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
  }
  
  Eigen::VectorXd get_support() const { return position_; }
  double get_support_min() const { return pos_min_; }
  double get_support_max() const { return pos_max_; }
  Eigen::SparseMatrix<double> get_binner() const { return binner_; }
  unsigned get_nbins() const { return nbins_; }
  Eigen::VectorXd get_nobs() const { return nobs_; }
  
  BINLESS_FORBID_COPY(SummarizerSettings);
  
private:
  Eigen::VectorXd position_;
  double pos_min_, pos_max_;
  Eigen::SparseMatrix<double> binner_; // Nbins x Ndata binary matrix
  unsigned nbins_;
  Eigen::VectorXd nobs_;
};

template<>
class FitterSettings<Bias,GAM> {
  
public:
  FitterSettings(const SummarizerSettings<BiasGAM>& settings, const BiasGAMConfig& conf) :
      max_iter_(conf.max_iter), tol_val_(conf.tol_val), sigma_(conf.sigma),
      K_(conf.bf_per_kb*(settings.get_support_max()-settings.get_support_min())/1000.),
      nbins_(settings.get_nbins()), nobs_(settings.get_nobs())  {
    auto position = settings.get_support();
    auto pos_min = settings.get_support_min();
    auto pos_max = settings.get_support_max();
    //X: design matrix
    X_ = generate_spline_base(position, pos_min, pos_max, get_K());
    //D: build difference matrix
    D_ = second_order_difference_matrix(get_K());
    //C: build constraint matrix to center
    if (conf.constraint_every > 0) {
      unsigned constraint_number = (pos_max-pos_min)/conf.constraint_every;
      Rcpp::Rcout << "Using " << constraint_number << "=" << (pos_max-pos_min) << "/" << conf.constraint_every << " constraints to model bias\n";
      const auto design2 = make_even_bins(pos_min, pos_max, constraint_number);
      const bool drop = true;
      Ceq_ = bin_data(position, design2, drop);
    } else {
      Ceq_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(0,K_);
    }
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
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, Cin); //not used
  BINLESS_GET_SET_DECL(Eigen::SparseMatrix<double>, const Eigen::SparseMatrix<double>&, Ceq);
  
  BINLESS_FORBID_COPY(FitterSettings);
};

template<>
struct FitterTraits<Bias,GAM> {
  typedef AnalyticalGAMLibrary library;
  static const bool has_inequality_constraints = false;
  static const bool has_equality_constraints = true;
};

//typedef Config<Bias,GAM> BiasConfig;
//typedef Estimator<Bias,GAM> BiasEstimator;

}
}

#endif

