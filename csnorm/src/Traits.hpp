#ifndef TRAITS_HPP
#define TRAITS_HPP

struct EstimatedOffset {};
struct ZeroOffset {};

struct PositiveSign {};
struct AnySign {};

struct BIC { enum { k = 0 }; };

template<unsigned kSD> struct CVkSD { enum { k = kSD }; };
typedef CVkSD<0> CV;


class SignalLikelihood;
class SignalRawData;
class SignalBinnedData;

struct Signal {
    typedef SignalLikelihood likelihood_t;
    typedef SignalRawData raw_t;
    typedef SignalBinnedData binned_t;
};

class DifferenceLikelihood;
class DifferenceRawData;
class DifferenceBinnedData;

struct Difference {
    typedef DifferenceLikelihood likelihood_t;
    typedef DifferenceRawData raw_t;
    typedef DifferenceBinnedData binned_t;
};

#endif

