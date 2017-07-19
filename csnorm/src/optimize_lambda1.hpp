#ifndef OPTIMIZE_LAMBDA1_HPP
#define OPTIMIZE_LAMBDA1_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

//objective functor to find lambda1 assuming eCprime=0, using BIC
struct obj_lambda1_BIC {
  obj_lambda1_BIC(double minUB, double tol_val,
                  IntegerVector patchno, NumericVector forbidden_vals,
                  NumericVector value, NumericVector weight, NumericVector valuehat,
                  NumericVector ncounts);
  
  double operator()(double x) const;
  
  NumericVector get(double val, std::string msg = "") const;
  
  double minUB_, minabsval_, maxabsval_, tol_val_, lsnc_;
  NumericVector forbidden_vals_;
  IntegerVector patchno_;
  NumericVector value_, absval_, weight_, valuehat_;
};

//objective functor to find lambda1 assuming eCprime=0, using CV
struct obj_lambda1_CV {
    obj_lambda1_CV(double minUB, double tol_val,
                IntegerVector patchno, NumericVector forbidden_vals,
                NumericVector value, NumericVector weight, NumericVector valuehat,
                NumericVector ncounts);

    double operator()(double x) const;

    NumericVector get(double val, std::string msg = "") const;

    double minUB_, minabsval_, maxabsval_, tol_val_, lsnc_;
    NumericVector forbidden_vals_;
    IntegerVector patchno_;
    NumericVector value_, absval_, weight_, valuehat_;
};

NumericVector cpp_optimize_lambda1(const DataFrame mat, int nbins,
                                   double tol_val, bool positive,
                                   double lambda1_min, int refine_num);


#endif

