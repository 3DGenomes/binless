#ifndef SCORE_PREPARATOR_HPP
#define SCORE_PREPARATOR_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "Traits.hpp"

template<typename Score, typename GaussianEstimator> class ScorePreparator;

// A class that produces data used for cross-validation of fused lasso regression results
template<int kSD, typename GaussianEstimator> class ScorePreparator<CVkSD<kSD>, GaussianEstimator> {
    
public:
    
    ScorePreparator(GaussianEstimator& gauss, const BinnedDataCore& binned, double lambda2, double y_default = -100.)
     : gauss_(gauss), binned_(binned), y_default_(y_default), lambda2_(lambda2), N_(binned.get_bin1().size()) {
        prepare();
        compute();
    }
    
    std::vector<int> get_cvgroup() const { return cvgroup_; }
    Rcpp::IntegerVector get_assembler_var() const { return Rcpp::wrap(cvgroup_); }
    
    std::vector<double> get_beta_cv() const { return beta_cv_; }
    Rcpp::NumericVector get_likelihood_var() const { return Rcpp::wrap(beta_cv_); }
    
    std::vector<double> get_betas() const { return betas_; }
    
private:
    
    void prepare();
    
    //run 2-fold cv calculation on checkerboard pattern
    void compute();
    
    GaussianEstimator& gauss_;
    const BinnedDataCore& binned_;
    const double y_default_, lambda2_;
    const unsigned N_;
    std::vector<int> cvgroup_;
    std::vector<double> beta_cv_;
    std::vector<std::vector<double> > betas_;
};

// A class that produces data used for BIC of fused lasso regression results
template<typename GaussianEstimator> class ScorePreparator<BIC, GaussianEstimator> {
    
public:
    
    ScorePreparator(GaussianEstimator& gauss, const BinnedDataCore& binned, double lambda2)
     : gauss_(gauss), binned_(binned), N_(binned.get_bin1().size())  {
        compute(binned.get_beta(), lambda2);
    }
    
    Rcpp::IntegerVector get_assembler_var() const { return binned_.get_ncounts(); }
    
    Rcpp::NumericVector get_likelihood_var() const { return Rcpp::wrap(beta_); }
   
private:
    
    //run 2-fold cv calculation on checkerboard pattern
    void compute(const Rcpp::NumericVector& beta_init, double lambda2);
    
    GaussianEstimator& gauss_;
    const BinnedDataCore& binned_;
    const unsigned N_;
    std::vector<double> beta_;
};

#include "ScorePreparator.ipp"

#endif

