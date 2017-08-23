#ifndef IRLS_ESTIMATOR_HPP
#define IRLS_ESTIMATOR_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp"

// A class that computes the 2D triangle grid fused lasso solution on some data.
// This class assumes that weights can change and uses an underlying gaussian model in an iterative way (IRLS).
template<typename GaussianEstimator, typename WeightsUpdater>
class IRLSEstimator {
    
public:
    
    struct var_t {
        GaussianEstimator::var_t gauss_var_;
        WeightsUpdater::var_t wt_var_;
        std::vector<double> get_beta() const { return gauss_var_.beta_; }
        set_beta(const std::vector<double>& beta) { gauss_var.beta_ = beta; }
    };
    
    //initialize the problem with a triangle grid with nrows
    //requesting precision to be below a given convergence criterion
    //final beta value will be clamped if clamp>0
    IRLSEstimator(unsigned nouter, double converge, GaussianEstimator& gauss, WeightsUpdater& wt)
    : nouter_(nouter), converge_(converge), counter_(0), gauss_(gauss), wt_(wt) {}
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // for a fixed lambda2. We alternate updating y and w using WeightsUpdater
    // and updating beta with the chosen fused lasso implementation in GaussianEstimator.
    // Initial state will be computed from beta_init
    void optimize(double lambda2, const var_t& init_data) {
        data_ = init_data;
        counter_ = 0;
        double precision = converge_+1;
        var_t data_old = data_;
        Rcpp::Rcout << " Perf iteration: start with lam2= " << lambda2 << " alpha= "
                    << gauss_.get_alpha() << " phi[0]= " << beta_[0] << "\n";
        do {
          //update weights
          wt_.update(data_.get_beta(), data_.wt_var_);
          auto y = wt_.get_y();
          auto w = wt_.get_w();
          //estimate beta
          gauss_.optimize(y, w, lambda2, data_.gauss_var_);
          data_.set_beta(gauss_.get());
          //update counters and compute precision
          precision = get_precision(data_.get_beta(),data_old.get_beta());
          ++counter_;
          data_old = data_;
          Rcpp::Rcout << " Iteration " << counter_ << " / " << nouter_ << " with lam2= " << lambda2 << " alpha= "
            << gauss_.get_alpha() << " reached maxval= " << precision
            << " after " << gauss_.get_ninner() << " steps " << " phi[0]= " << beta_[0] << "\n";
        } while (counter_ <= nouter_ && precision > converge_ );
        Rcpp::Rcout << " Perf iteration: end with lam2= " << lambda2 << " alpha= "
        << gauss_.get_alpha() << " phi[0]= " << beta_[0] << "\n";
        
    }
    
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
    
    GaussianEstimator& gauss_;
    WeightsUpdater& wt_;
    
    var_t data_;
    
};

//named constructor
template<typename GaussianEstimator, typename WeightsUpdater>
IRLSEstimator<GaussianEstimator, WeightsUpdater>
make_IRLSEstimator(unsigned nouter, double converge, GaussianEstimator& gauss, WeightsUpdater& wt) {
    return IRLSEstimator<GaussianEstimator, WeightsUpdater>(nouter,converge,gauss,wt);
}


#endif

