#ifndef OFFSET_HPP
#define OFFSET_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <utility> //pair

#include "Data.hpp"

struct EstimatedOffset {
  typedef std::pair<double,double> bounds_t;
  
  EstimatedOffset(const Data& data, double minval) :
       value_(data.get_value()), weight_(data.get_weight()), valuehat_(data.get_valuehat()), minval_(minval) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, weight_, valuehat_;
    double minval_;
};

struct ZeroOffset {
    typedef double bounds_t;
    
    ZeroOffset(const Data& data, double minUB) :
      value_(data.get_value()), absval_(abs(data.get_value())), weight_(data.get_weight()),
      valuehat_(data.get_valuehat()), minabsval_(min(absval_)), maxabsval_(max(absval_)), minUB_(minUB) {}

    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double UB) const;
    
private:
    NumericVector value_, absval_, weight_, valuehat_;
    double minabsval_, maxabsval_, minUB_;
};

#endif

