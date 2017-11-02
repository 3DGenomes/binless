#include <Rcpp.h>
#include <vector>

#include "util.hpp"
#include "Traits.hpp"

template<unsigned kSD, typename GaussianEstimator>
void ScorePreparator<CVkSD<kSD>,GaussianEstimator>::prepare() {
    //build cv groups
    const Rcpp::IntegerVector bin1(binned_.get_bin1()), bin2(binned_.get_bin2());
    const int ngroups=2;
    for (int i=0; i<N_; ++i)
        cvgroup_.push_back( (bin2(i)+bin1(i)) % ngroups ); // 2 cv groups in checkerboard pattern
    //resize beta_cv_
    beta_cv_ = std::vector<double>(N_,beta_default_);
    betas_.reserve(2);
}

template<unsigned kSD, typename GaussianEstimator>
void ScorePreparator<CVkSD<kSD>,GaussianEstimator>::compute() {
    //Compute fused lasso solutions on each group and report to beta_cv
    std::vector<double> y = Rcpp::as<std::vector<double> >(binned_.get_betahat());
    std::vector<double> w = Rcpp::as<std::vector<double> >(binned_.get_weight());
    const int ngroups=2;
    auto state = gauss_.get_state();
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
        gauss_.set_state(state);
        gauss_.optimize(y_cv, w_cv, lambda2_);
        std::vector<double> values = gauss_.get();
        betas_.push_back(values);
        
        //store fused solution at group positions back in beta_cv
        for (int i=0; i<N_; ++i) if (cvgroup_[i]==g) beta_cv_[i] = values[i];
    }
}
    

   
template<typename GaussianEstimator>
void ScorePreparator<BIC,GaussianEstimator>::compute(const Rcpp::NumericVector& beta_init, double lambda2) {
    if (beta_init.size() != N_) Rcpp::Rcout << "ERROR: wrong size for input to BIC calculation, check code\n";
    
    //Compute fused lasso solutions on each group and report to beta_cv
    std::vector<double> beta_init_v = Rcpp::as<std::vector<double> >(beta_init);
    std::vector<double> y = Rcpp::as<std::vector<double> >(binned_.get_betahat());
    std::vector<double> w = Rcpp::as<std::vector<double> >(binned_.get_weight());
    //compute fused lasso
    gauss_.optimize(y, beta_init_v, w, lambda2);
    beta_ = gauss_.get();
}

