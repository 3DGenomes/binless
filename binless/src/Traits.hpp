#ifndef TRAITS_HPP
#define TRAITS_HPP

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

#endif

