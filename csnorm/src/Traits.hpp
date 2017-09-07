#ifndef TRAITS_HPP
#define TRAITS_HPP

struct EstimatedOffset {};
struct ZeroOffset {};

struct PositiveSign {};
struct AnySign {};

struct BIC { enum { k = 0 }; };

template<unsigned kSD> struct CVkSD { enum { k = kSD }; };
typedef CVkSD<0> CV;


class SignalRawData;

struct Signal {
    typedef SignalRawData raw_t;
};

class DifferenceRawData;

struct Difference {
    typedef DifferenceRawData raw_t;
};

struct AllowDegeneracy {};
struct ForbidDegeneracy {};

#endif

