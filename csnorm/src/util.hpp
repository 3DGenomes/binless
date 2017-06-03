#ifndef UTIL_HPP
#define UTIL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#define SQUARE(x) ((x)*(x))

std::vector<double> soft_threshold(const std::vector<double>& beta, double eCprime, double lam1);

#endif

