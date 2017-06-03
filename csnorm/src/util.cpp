#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>

std::vector<double> soft_threshold(const std::vector<double>& beta, double eCprime, double lam1) {
  std::vector<double> phi;
  phi.reserve(beta.size());
  for (std::vector<double>::const_iterator it = beta.begin(); it != beta.end(); ++it) {
    double val = *it - eCprime;
    phi.push_back((val > 0) ? std::max(0., val - lam1) : std::min(0., val + lam1));
  }
  return phi;
}

