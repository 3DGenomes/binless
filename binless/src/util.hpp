#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
#include <vector>
#include <map>
#include <Eigen/Sparse>

#include "typedefs.hpp"
#include "BinnedData.hpp"

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
SpMat bin_data(const Eigen::VectorXd& data, const BinDesign& begin2bin, bool drop = true);

//Bin data evenly across Nbins bins and return in Nbins x Ndata sparse matrix form.
//by default, drops unused bins (rows with only zeroes)
SpMat bin_data_evenly(const Eigen::VectorXd& data, unsigned Nbins, bool drop = true);

//build dense triangular view of a collection of hi-c matrices
//arguments are taken as ordered factors whose length must correspond to the number of levels
//returned matrix is of size n_names * n_bins * (n_bins+1) / 2
Rcpp::DataFrame create_empty_matrix_cpp(Rcpp::IntegerVector names, Rcpp::IntegerVector bins);

#endif

