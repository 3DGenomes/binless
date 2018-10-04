#ifndef GFL_LIBRARY_TRIANGLE_IMPL_HPP
#define GFL_LIBRARY_TRIANGLE_IMPL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "GFLLibraryImpl.hpp"

std::vector<std::vector<int> > triangle_grid_chain(int nrows);

//policy class that implements fused lasso solution using the GFL library
class GFLLibraryTriangleImpl : public GFLLibraryImpl {
public:
    
    //init by cold start
    GFLLibraryTriangleImpl(int nrows) : GFLLibraryImpl(nrows, nrows), N_(nrows*(nrows+1)/2) {
        store_trails(nrows);
        reset();
    }
    
    //cold start
    /*virtual*/ void reset();
    
    //return internal state
    /*virtual*/ GFLState_t get_state() const;
    
    //warm start
    /*virtual*/ void set_state(const GFLState_t& state);
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    /*virtual*/ void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge);
    
private:
    
    void store_trails(int nrows);
    
    unsigned N_; //size of the fused lasso problem
    int ntrails_;
    std::vector<unsigned> trails_, breakpoints_;
    unsigned tsz_;
    
    std::vector<double> beta_impl_, z_, u_;
};

#endif

