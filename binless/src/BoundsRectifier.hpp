#ifndef BOUNDS_RECTIFIER_HPP
#define BOUNDS_RECTIFIER_HPP

#include <Rcpp.h>
#include <type_traits>

#include "Traits.hpp"
#include "BinnedData.hpp"

//the type of bounds when passed around
typedef std::pair<double,double> bounds_t;

//BoundsRectifier does tag dispatching by the type of Sign and Offset
template<typename Offset, typename Sign> class BoundsRectifier {};

//specialization in the case where the offset is estimated
template<typename Sign> class BoundsRectifier<EstimatedOffset, Sign> {
    static_assert(std::is_same<Sign, PositiveSign>::value,
                  "When offset is estimated, sign must be constrained to be positive!");
public:
    BoundsRectifier(const BinnedDataCore& data) : minbeta_(min(data.get_beta())) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t rectify_bounds(double lambda1, double minUB) const {
        double LB = minbeta_;
        double UB = LB + 2*lambda1;
        if (minUB > UB) UB = minUB;
        return bounds_t{LB,UB};
    }
    
private:
    const double minbeta_;
};


//specialization in the case where the offset is fixed
template<typename Sign>
class BoundsRectifier<ZeroOffset, Sign> {
public:
    BoundsRectifier(const BinnedDataCore&) {}
    
    //given an UB candidate (and implicit dof), find the adequate UB and LB
    bounds_t rectify_bounds(double lambda1, double minUB) const {
        double UB = lambda1;
        if (minUB > lambda1) UB = minUB;
        return bounds_t{-UB,UB};
    }
};

#endif

