#ifndef BASE_OBJECTIVES_HPP
#define BASE_OBJECTIVES_HPP

#include <Rcpp.h>
using namespace Rcpp;

#include "ScoreComputer.hpp"
#include "DataLikelihoods.hpp"

struct obj_lambda1_eCprime_base {
    typedef std::pair<double, double> bounds_t;
    
    obj_lambda1_eCprime_base(NumericVector value, NumericVector weight, NumericVector valuehat, double minval) :
       value_(value), weight_(weight), valuehat_(valuehat), minval_(minval) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, weight_, valuehat_;
    double minval_;
};

struct obj_lambda1_base {
    typedef double bounds_t;
    obj_lambda1_base(NumericVector value, NumericVector weight, NumericVector valuehat, double minUB) :
      value_(value), absval_(abs(value)), weight_(weight), valuehat_(valuehat), minabsval_(min(absval_)),
      maxabsval_(max(absval_)), minUB_(minUB) {}

    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, absval_, weight_, valuehat_;
    double minabsval_, maxabsval_, minUB_;
};




#endif

