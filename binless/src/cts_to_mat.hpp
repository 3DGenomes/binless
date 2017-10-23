#ifndef CTS_TO_MAT_HPP
#define CTS_TO_MAT_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Traits.hpp"

template<typename> class RawData;
template<typename> class BinnedData;

void remove_outliers(const std::vector<int>& bin1, const std::vector<int>& bin2,
                     std::vector<double>& phihat_var, List outliers);

void cts_to_signal_mat(const RawData<Signal>& raw, double eCprime, const Rcpp::NumericVector& beta_phi, BinnedData<Signal>& binned);

void cts_to_diff_mat(const RawData<Difference>& raw, const Rcpp::NumericVector& phi_ref, const Rcpp::NumericVector& beta_delta, BinnedData<Difference>& binned);


#endif
