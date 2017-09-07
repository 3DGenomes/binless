#ifndef PREPARATION_HPP
#define PREPARATION_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "Traits.hpp"

template<typename Calculation, typename GaussianEstimator> class Preparation;

// A class that produces data used for cross-validation of fused lasso regression results
template<int kSD, typename GaussianEstimator> class Preparation<CVkSD<kSD>, GaussianEstimator> {
    
public:
    
    Preparation(GaussianEstimator& gauss, const BinnedDataCore& binned, double y_default = -100.);
    
    //run 2-fold cv calculation on checkerboard pattern
    void compute(const std::vector<double>& beta_init, double lambda2);
    
    std::vector<int> get_cvgroup() const { return cvgroup_; }
    std::vector<int> get_assembler_var() const { return cvgroup_; }
    
    std::vector<double> get_beta_cv() const { return beta_cv_; }
    std::vector<double> get_likelihood_var() const { return beta_cv_; }
    
    std::vector<double> get_betas() const { return betas_; }
    
private:
    GaussianEstimator& gauss_;
    const BinnedDataCore& binned_;
    const double y_default_;
    unsigned N_;
    std::vector<int> cvgroup_;
    std::vector<double> beta_cv_;
    std::vector<std::vector<double> > betas_;
};

//named constructor
template<typename GaussianEstimator>
Preparation<CV,GaussianEstimator>
make_CVEstimator(GaussianEstimator& gauss, BinnedDataCore& binned, double y_default = -100) {
    return Preparation<CV,GaussianEstimator>(gauss,binned,y_default);
}

// A class that produces data used for BIC of fused lasso regression results
template<typename GaussianEstimator> class Preparation<BIC, GaussianEstimator> {
    
public:
    
    Preparation(GaussianEstimator&, const BinnedDataCore& binned) : binned_(binned) {}
    
    //run 2-fold cv calculation on checkerboard pattern
    void compute(const std::vector<double>& beta_init, double lambda2);
    
    std::vector<int> get_assembler_var() const { return binned_.get_ncounts(); }
    
    std::vector<double> get_likelihood_var() const { return beta_; }
   
private:
    GaussianEstimator& gauss_;
    const BinnedDataCore& binned_;
    std::vector<double> beta_;
};

//named constructor
template<typename GaussianEstimator>
Preparation<BIC,GaussianEstimator>
make_BICEstimator(GaussianEstimator& gauss, BinnedDataCore& binned) {
    return Preparation<BIC,GaussianEstimator>(gauss,binned);
}

#include "Preparation.ipp"

#endif

