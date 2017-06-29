#ifndef OPTIMIZE_LAMBDA1_ECPRIME_HPP
#define OPTIMIZE_LAMBDA1_ECPRIME_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

//objective functor to find lambda1 and eCprime assuming the signal is positive, using BIC
struct obj_lambda1_eCprime_BIC {
  obj_lambda1_eCprime_BIC(double minval, double maxval, double tol_val,
                      bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                      NumericVector value, NumericVector weight, NumericVector valuehat,
                      NumericVector ncounts, double lambda2);
  
  double operator()(double x) const;
  
  NumericVector get(double val, std::string msg = "") const;
  
  double minval_, valrange_, tol_val_, lsnc_, lambda2_;
  bool constrained_;
  IntegerVector patchno_;
  NumericVector forbidden_vals_, value_, weight_, valuehat_;
};

//objective functor to find lambda1 and eCprime assuming the signal is positive, using CV
struct obj_lambda1_eCprime_CV {
    obj_lambda1_eCprime_CV(double minval, double maxval, double tol_val,
                        bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                        NumericVector value, NumericVector weight, NumericVector valuehat,
                        NumericVector ncounts, double lambda2,
                        NumericVector beta_cv);

    double operator()(double x) const;

    NumericVector get(double val, std::string msg = "") const;

    double minval_, valrange_, tol_val_, lsnc_, lambda2_;
    bool constrained_;
    IntegerVector patchno_;
    NumericVector forbidden_vals_, value_, weight_, valuehat_;
};

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num, double lambda2, NumericVector beta_cv);


#endif

