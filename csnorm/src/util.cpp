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

//return ascending list of the values of all detected patches in the matrix
NumericVector get_patch_values(NumericVector value, IntegerVector patchno) {
    int npatches = max(patchno)+1;
    NumericVector unique_values(npatches); //patchno starts at 0
    for (int i=0; i<patchno.size(); ++i) unique_values(patchno(i)) = value(i);
    std::sort(unique_values.begin(), unique_values.end());
    return unique_values;
}

std::vector<double> compute_phi_ref(const std::vector<double>& delta_r,
                                    const std::vector<double>& phihat,
                                    const std::vector<double>& phihat_var, const std::vector<double>& phihat_ref,
                                    const std::vector<double>& phihat_var_ref) {
  const int N = delta_r.size();
  std::vector<double> phi_ref_r;
  phi_ref_r.reserve(N);
  for (int i=0; i<N; ++i) {
    double val;
    if (phihat_var_ref[i]==INFINITY && phihat_var[i]==INFINITY) {
      val=(phihat_ref[i]+phihat[i])/2;
    } else {
      val=(phihat_ref[i]/phihat_var_ref[i] + (phihat[i]-delta_r[i])/phihat_var[i])
      /(1/phihat_var_ref[i]+1/phihat_var[i]);
    }
    val = std::max(val,0.);
    phi_ref_r.push_back(val);
  }
  return phi_ref_r;
}


