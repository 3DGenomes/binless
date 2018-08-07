#ifndef GFL_LIBRARY_TRAPEZOIDAL_HPP
#define GFL_LIBRARY_TRAPEZOIDAL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Settings.hpp"

std::vector<std::vector<int> > trapezoidal_grid_chain(int nrows, int maxdiag);

//policy class that implements fused lasso solution using the GFL library
class GFLLibrary_trapezoidal {
public:
    
    typedef Rcpp::List GFLState_t;
    
    //init by cold start
    GFLLibrary_trapezoidal(unsigned nrows, unsigned maxdiag) : nrows_(nrows), maxdiag_(std::min(maxdiag,nrows)),
     N_(maxdiag_*(nrows_-maxdiag_+1) + (maxdiag_-1)*maxdiag_/2), counter_(0),
     inflate_(Settings<GFLLibrary_trapezoidal>::get_inflate()),
     ninner_(Settings<GFLLibrary_trapezoidal>::get_ninner()), alpha_(Settings<GFLLibrary_trapezoidal>::get_alpha()) {
        store_trails(nrows_, maxdiag_);
        reset();
    }
    
    //cold start
    void reset();
    
    //warm start
    void set_state(const GFLState_t& state);
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge);
    
    //return modified quantities
    std::vector<double> get_beta() const {
        return beta_tri_;
    }
    
    unsigned get_ninner() const {
        return counter_;
    }
    
    double get_alpha() const {
        return alpha_;
    }
    
    //return internal state
    GFLState_t get_state() const;
    
//protected:
    //to avoid direct destruction by user
    ~GFLLibrary_trapezoidal() {}
    
    std::vector<double> extract_trapeze(const std::vector<double>& vec) const;
    std::vector<double> fill_triangle(const std::vector<double>& y, const std::vector<double>& w) const;
    friend void test_trap(unsigned nrows, unsigned maxdiag);
    
private:
    
    void store_trails(int nrows, int maxdiag);
  
    
    unsigned nrows_, maxdiag_;
    unsigned N_; //size of the fused lasso problem
    unsigned counter_;
    double inflate_;
    int ninner_;
    int ntrails_;
    std::vector<int> trails_, breakpoints_;
    unsigned tsz_;
    
    double alpha_;
    std::vector<double> beta_tri_, beta_trap_, z_, u_;
};

void test_trap(unsigned nrows, unsigned maxdiag);



#endif

