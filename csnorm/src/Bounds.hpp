#ifndef BOUNDS_HPP
#define BOUNDS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <type_traits>

#include "ScoreComputer.hpp"
#include "DataLikelihoods.hpp"

struct PositiveSign {
    bool is_valid() const;
};

struct AnySign {
    bool is_valid() const;
};

struct EstimatedOffset {
  typedef std::pair<double,double> bounds_t;
  
  EstimatedOffset(const SignalData& data, double minval) :
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


//abstract class template
template<typename Offset, typename Sign>
struct BoundsComputer : public Offset {
    typedef typename Offset::bounds_t bounds_t;
    BoundsComputer(const Data& data, double bounds_specific);
};

//specialization in the case where the offset is estimated
template<typename Sign>
struct BoundsComputer<EstimatedOffset, Sign> : public EstimatedOffset {
    static_assert(std::is_same<Sign, PositiveSign>::value, "When offset is estimated, sign must be constrained to be positive!");
    BoundsComputer(const Data& data, double bounds_specific) : EstimatedOffset(data,bounds_specific) {}
};


//specialization in the case where the offset is fixed
template<typename Sign>
struct BoundsComputer<ZeroOffset, Sign> : public ZeroOffset {
    BoundsComputer(const Data& data, double bounds_specific) : ZeroOffset(data,bounds_specific) {}
};


#endif

