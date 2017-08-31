#ifndef BOUNDS_COMPUTER_HPP
#define BOUNDS_COMPUTER_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <type_traits>

#include "Offset.hpp"
#include "Sign.hpp"

class Data;

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

