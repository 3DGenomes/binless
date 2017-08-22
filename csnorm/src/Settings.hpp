#ifndef SETTINGS_HPP
#define SETTINGS_HPP

//class to store all settings used in the c++ side 
//to add settings for a new class, forward-declare the class
//and then specialize the Settings template on it
template<typename T> struct Settings {};

//settings for GFLLibrary
class GFLLibrary;
template<> struct Settings<GFLLibrary> {
    //inflation factor for adamts step update
    static const double get_inflate() { return 2.; }
    //maximum number of adamts steps
    static const int get_ninner() { return 100000; }
};

//settings for FusedLassoGaussianEstimator
template<typename Library> class FusedLassoGaussianEstimator;
template<typename Library> struct Settings<FusedLassoGaussianEstimator<Library> > {
    //at which value to clamp computed beta. Negative to turn off
    static const double get_clamp() { return 50; }
};

//settings for IRLSEstimator
template<typename Estimator, typename WeightsUpdater>
class IRLSEstimator;
template<typename Estimator, typename WeightsUpdater>
struct Settings<IRLSEstimator<Estimator, WeightsUpdater> > {
    //maximum number of outer IRLS iterations (weight updates)
    static const unsigned get_nouter() { return 500; }
};

#endif

