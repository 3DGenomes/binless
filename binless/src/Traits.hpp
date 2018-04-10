#ifndef TRAITS_HPP
#define TRAITS_HPP

/// TAGS

//Score
struct BIC { enum { k = 0 }; };
template<unsigned kSD> struct CVkSD { enum { k = kSD }; };
typedef CVkSD<0> CV;

//Offset
struct EstimatedOffset {};
struct ZeroOffset {};

//Sign
struct PositiveSign {};
struct AnySign {};

//Degeneracy
struct AllowDegeneracy {};
struct ForbidDegeneracy {};

//Calculation
struct Signal {};
struct Difference {};

namespace binless {
namespace fast {

//Leg
struct Decay {};
struct Bias {};

//Method
struct GAM {};
struct Mean {};


/// TRAITS declarations

template<typename Leg, typename Method>
struct Config;

template<typename Leg, typename Method>
struct FitterTraits;

template<typename Leg, typename Method>
struct FitterImpl;

class SummarizerSettings;

template<typename Leg, typename Method>
class SummarizerSettingsImpl;

template<typename Leg, typename Method>
class FitterSettings;

template<typename Method>
class Params;

}
}


#endif

