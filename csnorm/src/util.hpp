#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
#include <vector>

#include "BinnedData.hpp"

#define SQUARE(x) ((x)*(x))

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1);

Rcpp::NumericVector compute_phi_ref(const BinnedData<Difference>& binned,
                                    const Rcpp::NumericVector& delta);


#endif

