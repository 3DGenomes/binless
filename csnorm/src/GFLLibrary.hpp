#ifndef GFL_LIBRARY_HPP
#define GFL_LIBRARY_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

//policy class that implements fused lasso solution using the GFL library
class GFLLibrary {
public:
    
    GFLLibrary(unsigned nrows) : N_(nrows*(nrows+1)/2) {
        store_trails(nrows);
    }
    
    //temporary fix until the whole refactoring is done, and this moves to the constructor
    void setUp(double alpha, double inflate, int ninner, double converge, double clamp);
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    // values are clamped at their extreme value, see clamp parameter above
    // NOTE: for now, vectors are passed by value, to force copying from a conversion from Rcpp::NumericVector
    void optimize(const std::vector<double>& y, const std::vector<double>& beta_init,
                  const std::vector<double>& w, double lambda2) const;
    
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
    int ntrails_;
    std::vector<int> trails_, breakpoints_;
    unsigned tsz_;
    double inflate_;
    int ninner_;
    double converge_;
    double clamp_;
    
    //cached values are allowed to change in const call optimize()
    mutable double alpha_;
    mutable unsigned counter_;
    mutable std::vector<double> beta_;
};

#endif

