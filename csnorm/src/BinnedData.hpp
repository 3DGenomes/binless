#ifndef BINNED_DATA_HPP
#define BINNED_DATA_HPP

#include <Rcpp.h>

//generic binned data: input is betahat and weight, estimate is beta etc.
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
    
    Rcpp::NumericVector get_betahat() const { return betahat_; }
    void set_betahat(const Rcpp::NumericVector& betahat) { betahat_ = betahat; }
    
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
    Rcpp::NumericVector beta_, weight_, betahat_, ncounts_;
    Rcpp::IntegerVector patchno_, diag_idx_, diag_grp_;
};

//signal binned data: input is phihat and weight, estimate is beta, phi etc.
class SignalBinnedData : public BinnedData {
public:

    Rcpp::NumericVector get_beta_phi() const { return get_beta(); }
    void set_beta_phi(const Rcpp::NumericVector& beta_phi) { set_beta(beta_phi); }
    
    Rcpp::NumericVector get_phihat() const { return get_betahat(); }
    void set_phihat(const Rcpp::NumericVector& phihat) { set_betahat(phihat); }

    Rcpp::NumericVector get_phihat_var() const { return 1/get_weight(); }
};

//signal binned data: input is phihat and weight, estimate is beta, phi etc.
class DifferenceBinnedData : public BinnedData {
public:
    Rcpp::NumericVector get_beta_delta() const { return get_beta(); }
    void set_beta_delta(const Rcpp::NumericVector& beta_delta) { set_beta(beta_delta); }
    
    Rcpp::NumericVector get_deltahat() const { return get_betahat(); }
    void set_deltahat(const Rcpp::NumericVector& deltahat) { set_betahat(deltahat); }
    
    Rcpp::NumericVector get_phihat() const { return get_phihat_ref() + get_deltahat(); }

    Rcpp::NumericVector get_phihat_var() const { return 1/get_weight(); }
    
    Rcpp::NumericVector get_weight_ref() const { return weight_ref_; }
    void set_weight_ref(const Rcpp::NumericVector& weight_ref) { weight_ref_ = weight_ref; }
    
    Rcpp::NumericVector get_phihat_var_ref() const { return 1/get_weight_ref(); }

    Rcpp::NumericVector get_phihat_ref() const { return phihat_ref_; }
    void set_phihat_ref(const Rcpp::NumericVector& phihat_ref) { phihat_ref_ = phihat_ref; }
    
    Rcpp::NumericVector get_phi_ref() const { return phi_ref_; }
    void set_phi_ref(const Rcpp::NumericVector& phi_ref) { phi_ref_ = phi_ref; }
    
private:
    Rcpp::NumericVector weight_ref_, phihat_ref_, phi_ref_;
};



#endif

