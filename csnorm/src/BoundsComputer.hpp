#ifndef BOUNDS_COMPUTER_HPP
#define BOUNDS_COMPUTER_HPP

#include <Rcpp.h>
#include <type_traits>

#include "Traits.hpp"
#include "Data.hpp"

//the type of bounds when passed around
typedef std::pair<double,double> bounds_t;

//BoundsComputer does tag dispatching by the type of Sign and Offset
template<typename Offset, typename Sign> class BoundsComputer {};

//specialization in the case where the offset is estimated
template<typename Sign> class BoundsComputer<EstimatedOffset, Sign> {
    static_assert(std::is_same<Sign, PositiveSign>::value,
                  "When offset is estimated, sign must be constrained to be positive!");
public:
    BoundsComputer(const Data& data, double minval) :
     value_(data.get_value()), weight_(data.get_weight()),
     valuehat_(data.get_valuehat()), minval_(minval) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double val) const;
    
private:
    Rcpp::NumericVector value_, weight_, valuehat_;
    double minval_;
};


//specialization in the case where the offset is fixed
template<typename Sign>
class BoundsComputer<ZeroOffset, Sign> {
public:
    BoundsComputer(const Data& data, double minUB) :
    value_(data.get_value()), absval_(abs(data.get_value())), weight_(data.get_weight()),
    valuehat_(data.get_valuehat()), minabsval_(min(absval_)), maxabsval_(max(absval_)), minUB_(minUB) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double val) const;
    
private:
    Rcpp::NumericVector value_, absval_, weight_, valuehat_;
    double minabsval_, maxabsval_, minUB_;
};


#include "BoundsComputer.ipp" //implementation

#endif

