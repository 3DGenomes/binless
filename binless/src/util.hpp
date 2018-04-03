#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
#include <vector>
#include <map>
#include <Eigen/Sparse>

#include "BinnedData.hpp"

#define SQUARE(x) ((x)*(x))

//define a private variable and a public getter / setter pair
//CType is the type of the variable as passed in the setter (e.g. const Type&)
#define BINLESS_GET_SET_DECL(Type, CType, Name)                  \
public:                                                          \
  Type get_##Name() const {                                      \
    return Name##_;                                              \
  };                                                             \
  void set_##Name(CType value) {                                 \
    Name##_ = value;                                             \
  }                                                              \
  private:                                                       \
    Type Name##_;

//define a const private variable and public getter
#define BINLESS_GET_CONST_DECL(Type, Name)                         \
  public:                                                          \
    Type get_##Name() const {                                      \
      return Name##_;                                              \
    };                                                             \
    private:                                                       \
      const Type Name##_;
    
//define a private variable and public getter that returns a reference to it
#define BINLESS_GET_REF_DECL(Type, Name)                     \
public:                                                        \
Type& get_##Name() { return Name##_; }                         \
const Type& get_##Name() const { return Name##_; }             \
private:                                                       \
  Type Name##_;

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1);

Rcpp::NumericVector compute_phi_ref(const BinnedData<Difference>& binned,
                                    const Rcpp::NumericVector& delta);


double get_maximum_admissible_weight(const Rcpp::NumericVector& w, double percentile);

struct BinDesign {
  unsigned Nbins;
  std::map<double,unsigned> map;
  std::vector<double> begins, ends;
};

//create design for evenly sized bins with labels from 0 to N. Form is begin -> bin_no.
//Contains final bounday (size is N+1).
BinDesign make_even_bins(double lower, double upper, unsigned Nbins);

//bin data given a boundary map in the form [a,b) [b,c) ... [y,z] with labels from 0 to N-1
//Return in Nbins x Ndata sparse matrix form.
//by default, drops unused bins (rows with only zeroes)
Eigen::SparseMatrix<double, Eigen::RowMajor> bin_data(const Eigen::VectorXd& data, const BinDesign& begin2bin, bool drop = true);

//Bin data evenly across Nbins bins and return in Nbins x Ndata sparse matrix form.
//by default, drops unused bins (rows with only zeroes)
Eigen::SparseMatrix<double, Eigen::RowMajor> bin_data_evenly(const Eigen::VectorXd& data, unsigned Nbins, bool drop = true);


#endif

