#ifndef FUSED_LASSO_GAUSSIAN_ESTIMATOR_HPP
#define FUSED_LASSO_GAUSSIAN_ESTIMATOR_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "util.hpp"
#include "Settings.hpp"

// A class that computes the 2D triangle grid fused lasso solution on some data.
// This class uses a gaussian model, hence assumes that weights are held constant.
// This class defines the full interface and implements the sparse part of the lasso,
// while the Library policy contains the dense implementation
template<typename Library>
class FusedLassoGaussianEstimator : private Library {
    
public:
    
    //abstract type for variables optimized by Library. All must contain beta
    //each implementation can have additional parameters which will be carried along
    struct var_t : public Library::var_t {
        std::vector<double> beta_;
        //implicit conversion possible if other_var has beta member, at least
        template<typename T> var_t(const T& other_var) : Library::var_t(other_var), beta_(other_var.beta_) {}
    };
    
    //The problem to be computed is on a triangle grid with nrows
    //requesting precision to be below a given convergence criterion
    //final beta value will be clamped if clamp>0
    FusedLassoGaussianEstimator(unsigned nrows, double converge) : Library(nrows, converge),
    clamp_(Settings<FusedLassoGaussianEstimator<Library> >::get_clamp()) {}
    
    //Call this method to set starting values directly
    void initialize(const var_t& init_data) {
        GFLLibrary::initialize(init_data.beta_, init_data); //beta is passed separately
                                                            //because init_data gets downcasted
    }
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2,
                   const var_t& init_data) {
        Library::initialize(init_data);
        Library::optimize(y, w, lambda2);
    }
    
    std::vector<double> get_beta() const { Library::get_beta(); }
    
    //return soft-thresholded value, corresponding to the problem
    // sum_i w_i(y_i-offset-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j| + lambda1 * sum_i |beta_i|
    //where the solution beta for a given lambda2 was already computed
    //values will be clamped if necessary
    //returns an empty vector if optimize has not been called
    std::vector<double> get(double offset=0., double lambda1=0.) const {
        return soft_threshold(clamp(get_beta()), offset, lambda1);
    }
    
private:
    std::vector<double> clamp(std::vector<double> beta) const {
        //clamp values at +- clamp_ if needed
        if (clamp_>0) {
            for (double& i : beta) {
                i = std::min(clamp_, std::max(-clamp_, i));
            }
        }
        return beta;
    }
    
    const double clamp_;
    
};


#endif

