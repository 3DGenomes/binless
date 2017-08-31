#ifndef CALCULATION_TRAITS_HPP
#define CALCULATION_TRAITS_HPP

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

