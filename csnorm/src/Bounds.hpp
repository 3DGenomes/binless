#ifndef BOUNDS_HPP
#define BOUNDS_HPP

#include <Rcpp.h>
using namespace Rcpp;

#include "ScoreComputer.hpp"
#include "DataLikelihoods.hpp"

struct PositiveEstimatedOffsetBounds {
    typedef std::pair<double, double> bounds_t;
    
    PositiveEstimatedOffsetBounds(const SignalData& data, double minval) :
       value_(data.get_value()), weight_(data.get_weight()), valuehat_(data.get_valuehat()), minval_(minval) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, weight_, valuehat_;
    double minval_;
};

struct AnySignZeroOffsetBounds {
    typedef double bounds_t;
    
    AnySignZeroOffsetBounds(const Data& data, double minUB) :
      value_(data.get_value()), absval_(abs(data.get_value())), weight_(data.get_weight()),
      valuehat_(data.get_valuehat()), minabsval_(min(absval_)), maxabsval_(max(absval_)), minUB_(minUB) {}

    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, absval_, weight_, valuehat_;
    double minabsval_, maxabsval_, minUB_;
};




#endif

