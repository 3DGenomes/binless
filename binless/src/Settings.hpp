#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include "Traits.hpp"

//class to store all settings used in the c++ side 
//to add settings for a new class, forward-declare the class
//and then specialize the Settings template on it
template<typename T> struct Settings {};

//settings for GFLLibraryImpl
class GFLLibraryImpl;
template<> struct Settings<GFLLibraryImpl> {
    //initial adamts step size
    static const double get_alpha() { return 5.; }
    //inflation factor for adamts step update
    static const double get_inflate() { return 2.; }
    //maximum number of adamts steps
    static const int get_ninner_max() { return 30000; }
};


//settings for FusedLassoGaussianEstimator
template<typename Library> class FusedLassoGaussianEstimator;
template<typename Library> struct Settings<FusedLassoGaussianEstimator<Library> > {
    //at which value to clamp computed beta. Negative to turn off
    static const double get_clamp() { return 50; }
};

//settings for CV calculation
template<typename Score, typename GaussianEstimator> class ScorePreparator {};
template<unsigned kSD, typename GaussianEstimator> struct Settings<ScorePreparator<CVkSD<kSD>,GaussianEstimator> > {
    //set y to this value when left out during cv calculation. Has no influence except when lam2==0
    static const double get_y_default() { return -100; }
    //set beta_cv to this value upon initialization. All values will be overwritten
    static const double get_beta_default() { return -100; }
};

//settings for Likelihood
template<typename Calculation> class Likelihood {};
template<typename Calculation> struct Settings<Likelihood<Calculation> > {
    //clamp weights above this percentile to their maximum value
    static const double get_weights_max_percentile() { return 100; }
};


#endif

