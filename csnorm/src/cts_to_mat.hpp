#ifndef CTS_TO_MAT_HPP
#define CTS_TO_MAT_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

void remove_outliers(const std::vector<int>& bin1, const std::vector<int>& bin2,
                     std::vector<double>& phihat_var, List outliers);

DataFrame cts_to_signal_mat(const DataFrame& cts, int nbins, double dispersion,
                            const std::vector<double>& phi,
                            double eCprime, List outliers);

DataFrame cts_to_diff_mat(const DataFrame& cts, const DataFrame ref, int nbins,
                          double dispersion,
                          const std::vector<double>& phi_ref, const std::vector<double>& delta, List outliers);


#endif

