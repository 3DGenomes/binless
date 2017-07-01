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
    IntegerVector patchno_;
    NumericVector forbidden_vals_, absval_, value_, weight_, valuehat_;
};

NumericVector refine_minimum(const obj_lambda1_BIC& obj, double lam1,
                             double lam1_min, int refine_num, NumericVector patchvals, bool positive);

NumericVector cpp_optimize_lambda1(const DataFrame mat, int nbins,
                                   double tol_val, bool positive,
                                   double lambda1_min, int refine_num);


#endif

