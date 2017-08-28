#ifndef OPTIMIZE_LAMBDA1_ECPRIME_HPP
#define OPTIMIZE_LAMBDA1_ECPRIME_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <utility> //pair

struct obj_lambda1_eCprime_base {
    typedef std::pair<double, double> bounds_t;
    obj_lambda1_eCprime_base(NumericVector value, NumericVector weight, NumericVector valuehat, double minval);
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, weight_, valuehat_;
    double minval_;
};

//objective functor to find lambda1 and eCprime assuming the signal is positive, using BIC
struct obj_lambda1_eCprime_BIC : private obj_lambda1_eCprime_base {
  obj_lambda1_eCprime_BIC(double tol_val,
                      bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                      NumericVector value, NumericVector weight, NumericVector valuehat,
                      NumericVector ncounts, double lambda2);
  
  double operator()(double x) const;
  
  NumericVector get(double val, std::string msg = "") const;
  
  double tol_val_, lsnc_, lambda2_;
  bool constrained_;
  NumericVector forbidden_vals_;
  IntegerVector patchno_;
  NumericVector value_, weight_, valuehat_;
  double minval_;
};

//objective functor to find lambda1 and eCprime assuming the signal is positive, using CV
struct obj_lambda1_eCprime_CV : private obj_lambda1_eCprime_base {
    obj_lambda1_eCprime_CV(double minval, double tol_val,
                        bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                        NumericVector value, NumericVector weight, NumericVector valuehat,
                        NumericVector ncounts, double lambda2, IntegerVector cv_grp);

    double operator()(double x) const;

    NumericVector get(double val, std::string msg = "") const;

    double minval_, tol_val_, lsnc_, lambda2_;
    bool constrained_;
    NumericVector forbidden_vals_;
    IntegerVector patchno_;
    NumericVector value_, weight_, valuehat_;
    IntegerVector cv_grp_;
};

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num, double lambda2);


#endif

