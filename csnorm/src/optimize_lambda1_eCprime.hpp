#ifndef OPTIMIZE_LAMBDA1_ECPRIME_HPP
#define OPTIMIZE_LAMBDA1_ECPRIME_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

//objective functor to find lambda1 and eCprime assuming the signal is positive
struct obj_lambda1_eCprime {
    obj_lambda1_eCprime(double minval, double maxval, double tol_val,
                        bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                        NumericVector value, NumericVector weight, NumericVector valuehat,
                        NumericVector ncounts);

    double operator()(double x) const;

    NumericVector get(double lambda1, std::string msg = "") const;

    double minval_, valrange_, tol_val_, lsnc_;
    bool constrained_;
    IntegerVector patchno_;
    NumericVector forbidden_vals_, value_, weight_, valuehat_;
};

NumericVector refine_minimum(const obj_lambda1_eCprime& obj, double lam1,
                             double lam1_min, double lam1_max, int refine_num, NumericVector patchvals);

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num);


#endif

