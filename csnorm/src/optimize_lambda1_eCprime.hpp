#ifndef OPTIMIZE_LAMBDA1_ECPRIME_HPP
#define OPTIMIZE_LAMBDA1_ECPRIME_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <string>

//objective functor to find lambda1 and eCprime
struct obj_lambda1_eCprime {
  obj_lambda1_eCprime(double minval, double maxval, double tol_val,
                      bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                      NumericVector value, NumericVector weight, NumericVector valuehat,
                      NumericVector ncounts);
  
  double operator()(double const& x) const;
  
  NumericVector get(double const& lambda1, const std::string& msg) const;
  
  double minval_, maxval_, valrange_, tol_val_, lsnc_;
  bool constrained_;
  IntegerVector patchno_;
  NumericVector forbidden_vals_, value_, weight_, valuehat_;
};

NumericVector get_patch_values(NumericVector value, IntegerVector patchno);

NumericVector get_minimum_diagonal_values(NumericVector value, IntegerVector diag_idx);

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins, double tol_val, bool constrained,
                                           double lambda1_min);


#endif

