#ifndef SETTINGS_HPP
#define SETTINGS_HPP

//class to store all settings used in the c++ side 
template<typename T> struct Settings {};

//settings for GFLLibrary
class GFLLibrary;
template<> struct Settings<GFLLibrary> {
    static const double get_inflate() { return 2.; }
    static const int get_ninner() { return 100000; }
};

//settings for FusedLassoGaussianEstimator
template<typename Library> class FusedLassoGaussianEstimator;
template<typename Library> struct Settings<FusedLassoGaussianEstimator<Library> > {
    static const double get_clamp() { return 50; }
};

#endif

