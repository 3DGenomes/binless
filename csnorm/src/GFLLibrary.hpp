#ifndef GFL_LIBRARY_HPP
#define GFL_LIBRARY_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Settings.hpp"

//policy class that implements fused lasso solution using the GFL library
class GFLLibrary {
public:
    
    GFLLibrary(unsigned nrows, double converge) : N_(nrows*(nrows+1)/2), counter_(0), converge_(converge),
     inflate_(Settings<GFLLibrary>::get_inflate()), ninner_(Settings<GFLLibrary>::get_ninner()) {
        store_trails(nrows);
    }
    
    //temporary fix until the whole refactoring is done, and this moves to the constructor
    void setUp(double alpha);
    
    //prepare the next optimization round
    void prepare(const std::vector<double>& beta_init);
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2);
    
    //return modified quantities
    std::vector<double> get_beta() const {
        return beta_;
    }
    
    unsigned get_ninner() const {
        return counter_;
    }
    
    //temporary, see above
    double get_alpha() const {
        return alpha_;
    }
    
    
protected:
    //to avoid direct destruction by user
    ~GFLLibrary() {}
    
private:
    
    std::vector<std::vector<int> > triangle_grid_chain(unsigned nrows) const;
    
    void store_trails(unsigned nrows);
    
    unsigned N_; //size of the fused lasso problem
    unsigned counter_;
    double converge_; //convergence to this precision or better
    double inflate_;
    int ninner_;
    int ntrails_;
    std::vector<int> trails_, breakpoints_;
    unsigned tsz_;
    
    double alpha_;
    std::vector<double> beta_, z_, u_;
};

#endif

