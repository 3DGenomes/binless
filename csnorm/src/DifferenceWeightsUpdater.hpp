#ifndef DIFFERENCE_WEIGHTS_UPDATER_HPP
#define DIFFERENCE_WEIGHTS_UPDATER_HPP

#include <Rcpp.h>
#include <vector>

#include "util.hpp" //compute_phi_ref
#include "cts_to_mat.hpp" //cts_to_diff_mat
#include "RawData.hpp"

//policy class that updates IRLS data and weights given counts and a new beta
//This class takes care of the case of difference estimation
class DifferenceWeightsUpdater {
public:
    
    DifferenceWeightsUpdater(DifferenceRawData& data) : data_(data) {}
    
    void setUp(const std::vector<double>& phi_ref, const std::vector<double>& beta) {
        //store first estimate of phi_ref and compute mat
        Rcpp::DataFrame mat = cts_to_diff_mat(data_.get_cts(), data_.get_ref(), data_.get_nbins(), data_.get_dispersion(), phi_ref, beta, data_.get_outliers());
        data_.set_beta(beta);
        data_.set_phi_ref(phi_ref);
        data_.set_mat(mat);
    }
    
    void update(const std::vector<double>& beta) {
        //update phi_ref based on new beta
        DataFrame mat = data_.get_mat();
        std::vector<double> phihat = Rcpp::as<std::vector<double> >(mat["phihat"]);
        std::vector<double> phihat_var = Rcpp::as<std::vector<double> >
        (mat["phihat.var"]);
        std::vector<double> phihat_ref = Rcpp::as<std::vector<double> >
        (mat["phihat.ref"]);
        std::vector<double> phihat_var_ref = Rcpp::as<std::vector<double> >
        (mat["phihat.var.ref"]);
        std::vector<double> phi_ref = compute_phi_ref(beta, phihat, phihat_var, phihat_ref, phihat_var_ref);
        //compute new matrix of weights
        mat = cts_to_diff_mat(data_.get_cts(), data_.get_ref(), data_.get_nbins(), data_.get_dispersion(), phi_ref, beta, data_.get_outliers());
        data_.set_beta(beta);
        data_.set_phi_ref(phi_ref);
        data_.set_mat(mat);
    }
    
    std::vector<double> get_y() const {
        std::vector<double> phihat = Rcpp::as<std::vector<double> >(data_.get_mat()["phihat"]);
        const unsigned N = phihat.size();
        std::vector<double> y_r;
        y_r.reserve(N);
        std::vector<double> phi_ref = data_.get_phi_ref();
        for (int i=0; i<N; ++i) {
            y_r.push_back(phihat[i]-phi_ref[i]);
        }
        return y_r;
    }
    
    std::vector<double> get_w() const {
        std::vector<double> phihat_var = Rcpp::as<std::vector<double> >
        (data_.get_mat()["phihat.var"]);
        const unsigned N = phihat_var.size();
        std::vector<double> w_r;
        w_r.reserve(N);
        for (int i=0; i<N; ++i) {
            w_r.push_back(1/phihat_var[i]);
        }
        return w_r;
    }
    
    std::vector<int> get_bin1() const { return Rcpp::as<std::vector<int> > (data_.get_mat()["bin1"]); }
    
    std::vector<int> get_bin2() const { return Rcpp::as<std::vector<int> > (data_.get_mat()["bin2"]); }
    
    std::vector<double> get_phi_ref() const { return data_.get_phi_ref(); }
    
    Rcpp::DataFrame get_mat() const { return data_.get_mat(); }
    
private:
    DifferenceRawData& data_;
};

#endif

