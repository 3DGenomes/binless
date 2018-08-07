#ifndef GFL_LIBRARY_HPP
#define GFL_LIBRARY_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <memory>

#include "Settings.hpp"
#include "GFLLibraryTriangleImpl.hpp"
#include "GFLLibraryTrapezeImpl.hpp"

//policy class that implements fused lasso solution using the GFL library
//stores a handle to its implementation, which can be either triangle or trapezoidal
class GFLLibrary {
public:
    
    typedef Rcpp::List GFLState_t;
    
    //init by cold start
    GFLLibrary(unsigned nrows, unsigned maxdiag) {
        if (maxdiag > (nrows+1)/2 || maxdiag <= 1) {
          strategy_ = std::make_shared<GFLLibraryTriangleImpl>(nrows);
        } else {
          strategy_ = std::make_shared<GFLLibraryTrapezeImpl>(nrows, maxdiag);
        }
    }

    //init by cold start (triangle grid)
    GFLLibrary(unsigned nrows) : GFLLibrary(nrows,nrows) {}

    //reset
    void reset() { strategy_->reset(); }
    
    //return internal state
    GFLState_t get_state() const { return strategy_->get_state(); }
    
    //warm start
    void set_state(const GFLState_t& state) { strategy_->set_state(state); }
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
      strategy_->optimize(y,w,lambda2,converge);
    }
    
    //return modified quantities
    std::vector<double> get_beta() const { return strategy_->get_beta(); }
    
    unsigned get_ninner() const { return strategy_->get_ninner(); }

    double get_alpha() const { return strategy_->get_alpha(); }
    
//protected:
    //to avoid direct destruction by user (follows policy-based class design)
    ~GFLLibrary() = default;
    
private:
    std::shared_ptr<GFLLibraryImpl> strategy_;
};

#endif

