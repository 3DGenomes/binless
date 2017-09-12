#ifndef BOUNDS_COMPUTER_HPP
#define BOUNDS_COMPUTER_HPP

#include <Rcpp.h>
#include <type_traits>

#include "Traits.hpp"
#include "BinnedData.hpp"

//the type of bounds when passed around
typedef std::pair<double,double> bounds_t;

//BoundsComputer does tag dispatching by the type of Sign and Offset
template<typename Offset, typename Sign> class BoundsComputer {};

//specialization in the case where the offset is estimated
template<typename Sign> class BoundsComputer<EstimatedOffset, Sign> {
    static_assert(std::is_same<Sign, PositiveSign>::value,
                  "When offset is estimated, sign must be constrained to be positive!");
public:
    BoundsComputer(const BinnedDataCore& data) :
     beta_(data.get_beta()), weight_(data.get_weight()),
     y_(data.get_betahat()), minval_(min(data.get_beta())) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double val) const;
    
private:
    const Rcpp::NumericVector beta_, weight_, y_;
    const double minval_;
};


//specialization in the case where the offset is fixed
template<typename Sign>
class BoundsComputer<ZeroOffset, Sign> {
public:
    BoundsComputer(const BinnedDataCore& data) :
    beta_(data.get_beta()), absval_(abs(beta_)), weight_(data.get_weight()),
    y_(data.get_betahat()), minabsval_(min(absval_)), maxabsval_(max(absval_)) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t optimize_bounds(double val) const;
    
private:
    Rcpp::NumericVector beta_, absval_, weight_, y_;
    double minabsval_, maxabsval_;
};


#include "BoundsComputer.ipp" //implementation

#endif

