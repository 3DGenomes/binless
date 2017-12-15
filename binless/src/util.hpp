#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
#include <vector>
#include <map>
#include <Eigen/Sparse>

#include "BinnedData.hpp"

#define SQUARE(x) ((x)*(x))

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1);

Rcpp::NumericVector compute_phi_ref(const BinnedData<Difference>& binned,
                                    const Rcpp::NumericVector& delta);


double get_maximum_admissible_weight(const Rcpp::NumericVector& w, double percentile);

//create mapping for evenly sized bins with labels from 0 to N. Form is begin -> bin_no.
//Contains final bounday (size is N+1).
std::map<double,unsigned> make_even_bins(double lower, double upper, unsigned Nbins);

//bin data given a boundary map in the form [a,b) [b,c) ... [y,z] with labels from 0 to N-1
Eigen::SparseMatrix<unsigned> bin_data(const Eigen::VectorXd& data, const std::map<double,unsigned>& begin2bin);

//Bin data evenly across Nbins bins and return in Nbins x Ndata sparse matrix form.
Eigen::SparseMatrix<unsigned> bin_data_evenly(const Eigen::VectorXd& data, unsigned Nbins);



#endif

