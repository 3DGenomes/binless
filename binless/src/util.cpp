#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include "util.hpp"

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1) {
    std::vector<double> phi;
    phi.reserve(beta.size());
    for (std::vector<double>::const_iterator it = beta.begin(); it != beta.end();
            ++it) {
        double val = *it - eCprime;
        phi.push_back((val > 0) ? std::max(0., val - lam1) : std::min(0., val + lam1));
    }
    return phi;
}

Rcpp::NumericVector compute_phi_ref(const BinnedData<Difference>& binned,
                                    const Rcpp::NumericVector& delta) {
  const int N = delta.size();
  std::vector<double> phihat = Rcpp::as<std::vector<double> >(binned.get_phihat());
  std::vector<double> phihat_var = Rcpp::as<std::vector<double> >(binned.get_phihat_var());
  std::vector<double> phihat_ref = Rcpp::as<std::vector<double> >(binned.get_phihat_ref());
  std::vector<double> phihat_var_ref = Rcpp::as<std::vector<double> >(binned.get_phihat_var_ref());
  std::vector<double> phi_ref_r;
  phi_ref_r.reserve(N);
  for (int i=0; i<N; ++i) {
    double val;
    if (phihat_var_ref[i]==INFINITY && phihat_var[i]==INFINITY) {
      val=(phihat_ref[i]+phihat[i])/2;
    } else {
      val=(phihat_ref[i]/phihat_var_ref[i] + (phihat[i]-delta[i])/phihat_var[i])
      /(1/phihat_var_ref[i]+1/phihat_var[i]);
    }
    val = std::max(val,0.);
    phi_ref_r.push_back(val);
  }
  return Rcpp::wrap(phi_ref_r);
}

double get_maximum_admissible_weight(const Rcpp::NumericVector& w, double percentile) {
    //get sorted values of weights
    std::vector<double> sw = Rcpp::as<std::vector<double> >(w);
    std::sort(sw.begin(),sw.end());
    //find first nonzero weight
    unsigned i;
    for (i=0; i<sw.size(); ++i) if (sw[i]>0) break;
    //find percentile of positive values
    i = unsigned(percentile*(sw.size()-1-i)/100.)+i;
    //return its value
    return sw[i];
}