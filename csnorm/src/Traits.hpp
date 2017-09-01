#ifndef TRAITS_HPP
#define TRAITS_HPP

struct EstimatedOffset {};
struct ZeroOffset {};

struct PositiveSign {};
struct AnySign {};

class BIC { static const unsigned k = 0; };

template<unsigned kSD> class CVkSD { static const unsigned k = kSD; };
typedef CVkSD<0> CV;


class SignalLikelihood;
class SignalData;

struct Signal {
    typedef SignalLikelihood likelihood_t;
    typedef SignalData data_t;
};

class DifferenceLikelihood;
class DifferenceData;

struct Difference {
    typedef DifferenceLikelihood likelihood_t;
    typedef DifferenceData data_t;
};

#endif

