#ifndef IRLS_ESTIMATOR_HPP
#define IRLS_ESTIMATOR_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp"

// A class that computes the 2D triangle grid fused lasso solution on some data.
// This class assumes that weights can change and uses an underlying gaussian model in an iterative way (IRLS).
template<typename GaussianEstimator, typename WeightsUpdater>
class IRLSEstimator : public GaussianEstimator, public WeightsUpdater {
    
public:
    
    //initialize the problem with a triangle grid with nrows
    //requesting precision to be below a given convergence criterion
    //final beta value will be clamped if clamp>0
    IRLSEstimator(unsigned nrows, double converge, unsigned nouter, double dispersion,
                  const Rcpp::DataFrame& cts, const Rcpp::List& outliers)
    : GaussianEstimator(nrows, converge), WeightsUpdater(nrows, dispersion, cts, outliers),
      nouter_(nouter), converge_(converge), counter_(0) {}
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // for a fixed lambda2. We alternate updating y and w using WeightsUpdater
    // and updating beta with the chosen fused lasso implementation in GaussianEstimator.
    // Initial state will be computed from beta_init
    void optimize(const std::vector<double>& beta_init, double lambda2) {
        beta_ = beta_init;
        counter_ = 0;
        double precision = converge_+1;
        std::vector<double> beta_old = beta_;
        /*Rcpp::Rcout << " Perf iteration: start with lam2= " << lambda2 << " alpha= "
                    << GaussianEstimator::get_alpha() << " phi[0]= " << beta_[0] << "\n";*/
        do {
          //update weights
          WeightsUpdater::update(beta_);
          auto y = WeightsUpdater::get_y();
          auto w = WeightsUpdater::get_w();
          //estimate beta
          GaussianEstimator::optimize(y, beta_, w, lambda2);
          beta_ = GaussianEstimator::get();
          //update counters and compute precision
          precision = get_precision(beta_,beta_old);
          ++counter_;
          beta_old = beta_;
          /*Rcpp::Rcout << " Iteration " << counter_ << " with lam2= " << lambda2 << " alpha= "
            << GaussianEstimator::get_alpha() << " reached maxval= " << precision
            << " after " << GaussianEstimator::get_ninner() << " steps " << " phi[0]= " << beta_[0] << "\n";*/
        } while (counter_ <= nouter_ && precision > converge_ );
        /*Rcpp::Rcout << " Perf iteration: end with lam2= " << lambda2 << " alpha= "
        << GaussianEstimator::get_alpha() << " phi[0]= " << beta_[0] << "\n";*/
        
    }
    
    //the get() call is already inherited from GaussianEstimator
    //std::vector<double> get(double offset=0., double lambda1=0.) const
    
    //return the number of outer iterations
    unsigned get_nouter() const {
        return counter_;
    }
    
private:
    double get_precision(const std::vector<double>& beta, const std::vector<double>& beta_old) const {
        double maxval = std::abs(beta[0]-beta_old[0]);
        const unsigned N = beta.size();
        for (unsigned i=1; i < N; ++i) {
            maxval = std::max(std::abs(beta[i]-beta_old[i]), maxval);
        }
        return maxval;
    }
        
    const unsigned nouter_;
    const double converge_;
    
    unsigned counter_;
    std::vector<double> beta_;
    
};


#endif

