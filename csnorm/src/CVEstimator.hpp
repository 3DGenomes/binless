#ifndef CV_ESTIMATOR_HPP
#define CV_ESTIMATOR_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp"

// A class that produces data used for cross-validation of fused lasso regression results
template<typename GaussianEstimator, typename WeightsUpdater>
class CVEstimator {
    
public:
    
    CVEstimator(GaussianEstimator& gauss, WeightsUpdater& wt, double y_default = -100.)
      : gauss_(gauss), wt_(wt), y_default_(y_default), N_(wt.get_bin1().size()) {
        //build cv groups
        const std::vector<int> bin1(wt.get_bin1()), bin2(wt.get_bin2());
        const int ngroups=2;
        for (int i=0; i<N_; ++i)
            cvgroup_.push_back( (bin2[i]+bin1[i]) % ngroups ); // 2 cv groups in checkerboard pattern
        //resize beta_cv_
        const double beta_default = -100.; //all values will be overwritten anyways
        beta_cv_ = std::vector<double>(N_,beta_default);
        betas_.reserve(2);
    }
    
    //run 2-fold cv calculation on checkerboard pattern
    void compute(const std::vector<double>& beta_init, double lambda2) {
        if (beta_init.size() != N_) Rcpp::Rcout << "ERROR: wrong size for input to CV calculation, check code\n";
        
        //Compute fused lasso solutions on each group and report to beta_cv
        std::vector<double> y = wt_.get_y();
        std::vector<double> w = wt_.get_w();
        const int ngroups=2;
        for (int g=0; g<ngroups; ++g) {
            //prepare data and weights for group g and copy initial values
            std::vector<double> y_cv(y), w_cv(w);
            for (int i=0; i<N_; ++i) {
                if (cvgroup_[i]==g) {
                    y_cv[i]=y_default_; //essential if lam2==0
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
    
    std::vector<int> get_cvgroup() const { return cvgroup_; }
    
    std::vector<double> get_beta_cv() const { return beta_cv_; }
    
    std::vector<double> get_betas() const { return betas_; }
    
private:
    GaussianEstimator& gauss_;
    WeightsUpdater& wt_;
    const double y_default_;
    unsigned N_;
    std::vector<int> cvgroup_;
    std::vector<double> beta_cv_;
    std::vector<std::vector<double> > betas_;
};

//named constructor
template<typename GaussianEstimator, typename WeightsUpdater>
CVEstimator<GaussianEstimator, WeightsUpdater>
make_CVEstimator(GaussianEstimator& gauss, WeightsUpdater& wt, double y_default = -100) {
    return CVEstimator<GaussianEstimator, WeightsUpdater>(gauss,wt,y_default);
}

#endif

