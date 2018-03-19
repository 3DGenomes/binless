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

template<typename Calculation>
void cts_to_mat(const RawData<Calculation>& raw, BinnedData<Calculation>& binned);

Rcpp::DataFrame rcpp_cts_to_signal_mat(int nbins, double alpha, Rcpp::DataFrame cts, Rcpp::List metadata);

#endif

