#ifndef TAGS_HPP
#define TAGS_HPP

struct EstimatedOffset {};
struct ZeroOffset {};

struct PositiveSign {};
struct AnySign {};

class BIC {
    static const unsigned k = 0;
};

template<unsigned kSD> class CVkSD {
    static const unsigned k = kSD;
};
typedef CVkSD<0> CV;


#endif

