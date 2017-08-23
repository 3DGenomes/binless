#ifndef DIFFERENCE_WEIGHTS_UPDATER_HPP
#define DIFFERENCE_WEIGHTS_UPDATER_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp" //compute_phi_ref
#include "cts_to_mat.hpp" //cts_to_diff_mat

//policy class that updates IRLS data and weights given counts and a new beta
//This class takes care of the case of difference estimation
class DifferenceWeightsUpdater {
public:
    
    DifferenceWeightsUpdater(unsigned nrows, double dispersion,
                         const Rcpp::DataFrame& cts, const Rcpp::DataFrame& ref, const Rcpp::List& outliers) :
       nrows_(nrows), dispersion_(dispersion), cts_(cts), ref_(ref), outliers_(outliers) {}
    
    void setUp(const std::vector<double>& phi_ref, const std::vector<double>& beta) {
        //store first estimate of phi_ref and compute mat
        phi_ref_ = phi_ref;
        mat_ = cts_to_diff_mat(cts_, ref_, nrows_, dispersion_, phi_ref_, beta, outliers_);
    }
    
    void update(const std::vector<double>& beta) {
        //update phi_ref based on new beta
        std::vector<double> phihat = Rcpp::as<std::vector<double> >(mat_["phihat"]);
        std::vector<double> phihat_var = Rcpp::as<std::vector<double> >
        (mat_["phihat.var"]);
        std::vector<double> phihat_ref = Rcpp::as<std::vector<double> >
        (mat_["phihat.ref"]);
        std::vector<double> phihat_var_ref = Rcpp::as<std::vector<double> >
        (mat_["phihat.var.ref"]);
        phi_ref_ = compute_phi_ref(beta, phihat, phihat_var, phihat_ref, phihat_var_ref);
        //compute new matrix of weights
        mat_ = cts_to_diff_mat(cts_, ref_, nrows_, dispersion_, phi_ref_, beta, outliers_);
     }
    
    std::vector<double> get_y() const {
        std::vector<double> phihat = Rcpp::as<std::vector<double> >(mat_["phihat"]);
        const unsigned N = phihat.size();
        std::vector<double> y_r;
        y_r.reserve(N);
        for (int i=0; i<N; ++i) {
            y_r.push_back(phihat[i]-phi_ref_[i]);
        }
        return y_r;
    }
    
    std::vector<double> get_w() const {
        std::vector<double> phihat_var = Rcpp::as<std::vector<double> >
        (mat_["phihat.var"]);
        const unsigned N = phihat_var.size();
        std::vector<double> w_r;
        w_r.reserve(N);
        for (int i=0; i<N; ++i) {
            w_r.push_back(1/phihat_var[i]);
        }
        return w_r;
    }
    
    std::vector<double> get_phi_ref() const { return phi_ref_; }
    Rcpp::DataFrame get_mat() const { return mat_; }
    
private:
    const unsigned nrows_;
    const double dispersion_;
    const Rcpp::DataFrame cts_, ref_;
    const Rcpp::List outliers_;
    
    std::vector<double> phi_ref_;
    Rcpp::DataFrame mat_;
    
};

#endif

