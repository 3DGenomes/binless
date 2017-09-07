#ifndef TRAITS_HPP
#define TRAITS_HPP

struct EstimatedOffset {};
struct ZeroOffset {};

struct PositiveSign {};
struct AnySign {};

struct BIC { enum { k = 0 }; };

template<unsigned kSD> struct CVkSD { enum { k = kSD }; };
typedef CVkSD<0> CV;

struct Signal {};

struct Difference {};

struct AllowDegeneracy {};
struct ForbidDegeneracy {};

#endif

