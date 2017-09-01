#ifndef BINNED_DATA_HPP
#define BINNED_DATA_HPP

#include <Rcpp.h>

//generic binned data: input is y and weight, estimate is beta etc.
class BinnedData {
public:
    Rcpp::IntegerVector get_bin1() const { return bin1_; }
    void set_bin1(const Rcpp::IntegerVector& bin1) { bin1_ = bin1; }
    
    Rcpp::IntegerVector get_bin2() const { return bin2_; }
    void set_bin2(const Rcpp::IntegerVector& bin2) { bin2_ = bin2; }
    
    Rcpp::NumericVector get_beta() const { return beta_; }
    void set_beta(const Rcpp::NumericVector& beta) { beta_ = beta; }
    
    Rcpp::NumericVector get_weight() const { return weight_; }
    void set_weight(const Rcpp::NumericVector& weight) { weight_ = weight; }
    
    Rcpp::NumericVector get_y() const { return y_; }
    void set_y(const Rcpp::NumericVector& y) { y_ = y; }
    
    Rcpp::NumericVector get_ncounts() const { return ncounts_; }
    void set_ncounts(const Rcpp::NumericVector& ncounts) { ncounts_ = ncounts; }
    
    Rcpp::IntegerVector get_patchno() const { return patchno_; }
    void set_patchno(const Rcpp::IntegerVector& patchno) { patchno_ = patchno; }
    
    Rcpp::IntegerVector get_diag_idx() const { return diag_idx_; }
    void set_diag_idx(const Rcpp::IntegerVector& diag_idx) { diag_idx_ = diag_idx; }
    
    Rcpp::IntegerVector get_diag_grp() const { return diag_grp_; }
    void set_diag_grp(const Rcpp::IntegerVector& diag_grp) { diag_grp_ = diag_grp; }
    
private:
    Rcpp::IntegerVector bin1_, bin2_;
    Rcpp::NumericVector beta_, weight_, y_, ncounts_;
    Rcpp::IntegerVector patchno_, diag_idx_, diag_grp_;
};

//signal binned data: input is phihat and weight, estimate is beta, phi etc.
class SignalBinnedData : private BinnedData {
public:
    using BinnedData::get_bin1;
    using BinnedData::set_bin1;
    
    using BinnedData::get_bin2;
    using BinnedData::set_bin2;
    
    Rcpp::NumericVector get_beta_phi() const { return get_beta(); }
    void set_beta_phi(const Rcpp::NumericVector& beta_phi) { set_beta(beta_phi); }
    
    using BinnedData::get_weight;
    using BinnedData::set_weight;
    
    Rcpp::NumericVector get_phihat() const { return get_y(); }
    void set_phihat(const Rcpp::NumericVector& phihat) { set_y(phihat); }
    
    using BinnedData::get_ncounts;
    using BinnedData::set_ncounts;
    
    using BinnedData::get_patchno;
    using BinnedData::set_patchno;
    
    using BinnedData::get_diag_idx;
    using BinnedData::set_diag_idx;
    
    using BinnedData::get_diag_grp;
    using BinnedData::set_diag_grp;
};

//signal binned data: input is phihat and weight, estimate is beta, phi etc.
class DifferenceBinnedData : private BinnedData {
public:
    using BinnedData::get_bin1;
    using BinnedData::set_bin1;
    
    using BinnedData::get_bin2;
    using BinnedData::set_bin2;
    
    Rcpp::NumericVector get_beta_delta() const { return get_beta(); }
    void set_beta_delta(const Rcpp::NumericVector& beta_delta) { set_beta(beta_delta); }
    
    using BinnedData::get_weight;
    using BinnedData::set_weight;
    
    Rcpp::NumericVector get_deltahat() const { return get_y(); }
    void set_deltahat(const Rcpp::NumericVector& deltahat) { set_y(deltahat); }
    
    Rcpp::NumericVector get_phihat() const { return get_phihat_ref() + get_deltahat(); }
    
    using BinnedData::get_ncounts;
    using BinnedData::set_ncounts;
    
    using BinnedData::get_patchno;
    using BinnedData::set_patchno;
    
    using BinnedData::get_diag_idx;
    using BinnedData::set_diag_idx;
    
    using BinnedData::get_diag_grp;
    using BinnedData::set_diag_grp;
    
    Rcpp::NumericVector get_weight_ref() const { return weight_ref_; }
    void set_weight_ref(const Rcpp::NumericVector& weight_ref) { weight_ref_ = weight_ref; }
    
    Rcpp::NumericVector get_phihat_ref() const { return phihat_ref_; }
    void set_phihat_ref(const Rcpp::NumericVector& phihat_ref) { phihat_ref_ = phihat_ref; }
    
    Rcpp::NumericVector get_phi_ref() const { return phi_ref_; }
    void set_phi_ref(const Rcpp::NumericVector& phi_ref) { phi_ref_ = phi_ref; }
    
private:
    Rcpp::NumericVector weight_ref_, phihat_ref_, phi_ref_;
};



#endif

