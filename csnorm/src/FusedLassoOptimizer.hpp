#ifndef FUSED_LASSO_OPTIMIZER_HPP
#define FUSED_LASSO_OPTIMIZER_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "util.hpp"

// A class that computes the 2D triangle grid fused lasso solution on some data.
// This class defines the full interface and implements the sparse part of the lasso,
// while the Library policy contains the dense implementation
template<typename Library>
class FusedLassoOptimizer : public Library {
    
public:
    
    //initialize the problem with a triangle grid with nrows
    FusedLassoOptimizer(unsigned nrows) : Library(nrows) {}
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& beta_init,
                  const std::vector<double>& w, double lambda2) {
        Library::optimize(y, beta_init, w, lambda2);
    }
    
    //return soft-thresholded value, corresponding to the problem
    // sum_i w_i(y_i-offset-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j| + lambda1 * sum_i |beta_i|
    //where the solution beta for a given lambda2 was already computed
    //returns an empty vector if optimize has not been called
    std::vector<double> get(double offset=0., double lambda1=0.) const {
        return soft_threshold(Library::get_beta(), offset, lambda1);
    }
    
};


#endif

