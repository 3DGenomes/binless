#ifndef OPTIMIZE_LAMBDA1_DIFF_HPP
#define OPTIMIZE_LAMBDA1_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <utility> //pair


struct obj_lambda1_diff_base {
    typedef double bounds_t;
    obj_lambda1_diff_base(NumericVector value, NumericVector weight, NumericVector valuehat, double minUB);
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, absval_, weight_, valuehat_;
    double minabsval_, maxabsval_, minUB_;
};

//objective functor to find lambda1 assuming eCprime=0, using BIC
struct obj_lambda1_diff_BIC : private obj_lambda1_diff_base {
  obj_lambda1_diff_BIC(double minUB, double tol_val,
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
struct obj_lambda1_diff_CV : private obj_lambda1_diff_base {
    obj_lambda1_diff_CV(double minUB, double tol_val,
                IntegerVector patchno, NumericVector forbidden_vals,
                NumericVector value, NumericVector weight, NumericVector valuehat,
                NumericVector weight_ref, NumericVector valuehat_ref,
                NumericVector ncounts, IntegerVector cv_grp);

    double operator()(double x) const;

    NumericVector get(double val, std::string msg = "") const;

    double minUB_, minabsval_, maxabsval_, tol_val_, lsnc_;
    NumericVector forbidden_vals_;
    IntegerVector patchno_;
    NumericVector value_, absval_, weight_, valuehat_, weight_ref_, valuehat_ref_;
    IntegerVector cv_grp_;
};

NumericVector cpp_optimize_lambda1_diff(const DataFrame mat, int nbins,
                                   double tol_val, double lambda1_min, int refine_num);


#endif

