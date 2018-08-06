#ifndef GFL_LIBRARY_TRIANGLE_HPP
#define GFL_LIBRARY_TRIANGLE_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Settings.hpp"

std::vector<std::vector<int> > triangle_grid_chain(int nrows);

//policy class that implements fused lasso solution using the GFL library
class GFLLibrary_triangle {
public:
    
    typedef Rcpp::List GFLState_t;
    
    //init by cold start
    GFLLibrary_triangle(unsigned nrows) : N_(nrows*(nrows+1)/2), counter_(0),
     inflate_(Settings<GFLLibrary_triangle>::get_inflate()), ninner_(Settings<GFLLibrary_triangle>::get_ninner()),
     alpha_(Settings<GFLLibrary_triangle>::get_alpha()) {
        store_trails(nrows);
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
        return beta_;
    }
    
    unsigned get_ninner() const {
        return counter_;
    }
    
    double get_alpha() const {
        return alpha_;
    }
    
    //return internal state
    GFLState_t get_state() const;
    
protected:
    //to avoid direct destruction by user
    ~GFLLibrary_triangle() {}
    
private:
    
    void store_trails(int nrows);
    
    unsigned N_; //size of the fused lasso problem
    unsigned counter_;
    double inflate_;
    int ninner_;
    int ntrails_;
    std::vector<int> trails_, breakpoints_;
    unsigned tsz_;
    
    double alpha_;
    std::vector<double> beta_, z_, u_;
};

#endif

