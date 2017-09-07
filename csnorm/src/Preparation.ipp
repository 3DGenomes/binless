#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "Traits.hpp"

template<int kSD, typename GaussianEstimator>
Preparation<kSD,GaussianEstimator>::Preparation(GaussianEstimator& gauss, const BinnedDataCore& binned,
                                   double y_default = -100.)
      : gauss_(gauss), binned_(binned), y_default_(y_default), N_(binned.get_bin1().size()) {
    //build cv groups
    const Rcpp::IntegerVector bin1(binned.get_bin1()), bin2(binned.get_bin2());
    const int ngroups=2;
    for (int i=0; i<N_; ++i)
        cvgroup_.push_back( (bin2(i)+bin1(i)) % ngroups ); // 2 cv groups in checkerboard pattern
    //resize beta_cv_
    const double beta_default = -100.; //all values will be overwritten anyways
    beta_cv_ = std::vector<double>(N_,beta_default);
    betas_.reserve(2);
}
    
template<int kSD, typename GaussianEstimator>
void Preparation<kSD,GaussianEstimator>::compute(const std::vector<double>& beta_init, double lambda2) {
    if (beta_init.size() != N_) Rcpp::Rcout << "ERROR: wrong size for input to CV calculation, check code\n";
    
    //Compute fused lasso solutions on each group and report to beta_cv
    std::vector<double> y = Rcpp::as<std::vector<double> >(binned_.get_betahat());
    std::vector<double> w = Rcpp::as<std::vector<double> >(binned_.get_weight());
    const int ngroups=2;
    for (int g=0; g<ngroups; ++g) {
        //prepare data and weights for group g and copy initial values
        std::vector<double> y_cv(y), w_cv(w);
        for (int i=0; i<N_; ++i) {
            if (cvgroup_[i]==g) {
                y_cv[i]=y_default_; //has influence if lam2==0
                w_cv[i]=0;
            }
        }
        //compute fused lasso
        gauss_.optimize(y_cv, beta_init, w_cv, lambda2);
        std::vector<double> values = gauss_.get();
        betas_.push_back(values);
        
        //store fused solution at group positions back in beta_cv
        for (int i=0; i<N_; ++i) if (cvgroup_[i]==g) beta_cv_[i] = values[i];
    }
}
    

   
template<typename GaussianEstimator>
void Preparation<BIC,GaussianEstimator>::compute(const std::vector<double>& beta_init, double lambda2) {
    if (beta_init.size() != N_) Rcpp::Rcout << "ERROR: wrong size for input to BIC calculation, check code\n";
    
    //Compute fused lasso solutions on each group and report to beta_cv
    std::vector<double> y = Rcpp::as<std::vector<double> >(binned_.get_betahat());
    std::vector<double> w = Rcpp::as<std::vector<double> >(binned_.get_weight());
    //compute fused lasso
    gauss_.optimize(y, beta_init, w, lambda2);
    beta_ = gauss_.get();
}

