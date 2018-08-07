#ifndef GFL_LIBRARY_TRAPEZE_IMPL_HPP
#define GFL_LIBRARY_TRAPEZE_IMPL_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "GFLLibraryImpl.hpp"

std::vector<std::vector<int> > trapezoidal_grid_chain(int nrows, int maxdiag);

//policy class that implements fused lasso solution using the GFL library
class GFLLibraryTrapezeImpl : public GFLLibraryImpl {
public:
    
    //init by cold start
    GFLLibraryTrapezeImpl(int nrows, int maxdiag) : GFLLibraryImpl(nrows, maxdiag),
        nrows_(nrows), maxdiag_(maxdiag), N_(maxdiag_*(nrows_-maxdiag_+1) + (maxdiag_-1)*maxdiag_/2) {
        store_trails(nrows, maxdiag);
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
    
    void store_trails(int nrows, int maxdiag);
  
    std::vector<double> extract_trapeze(const std::vector<double>& vec) const;
    std::vector<double> fill_triangle(const std::vector<double>& y, const std::vector<double>& w) const;
  
    int nrows_, maxdiag_;
    unsigned N_; //size of the fused lasso problem
    int ntrails_;
    std::vector<int> trails_, breakpoints_;
    unsigned tsz_;
    
    std::vector<double> beta_impl_, z_, u_;
};

#endif

