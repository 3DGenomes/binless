#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <set>

#include "cts_to_mat.hpp"

#include "cts_core.h" //cts_to_signal_mat_core

#include "RawData.hpp"
#include "BinnedData.hpp"


void remove_outliers(const std::vector<int>& bin1, const std::vector<int>& bin2, std::vector<double>& phihat_var, List outliers) {
  unsigned nbetas = phihat_var.size();
  //remove bad rows
  const std::vector<int> bad_rows = as<std::vector<int> >(outliers["bad.rows"]);
  const std::set<int> bad_rows_set(bad_rows.begin(), bad_rows.end());
  for (unsigned i=0; i<nbetas; ++i) {
    if ( (bad_rows_set.find(bin1[i]) != bad_rows_set.end()) ||
         (bin2[i]!=bin1[i] && (bad_rows_set.find(bin2[i]) != bad_rows_set.end()) )
       ) { phihat_var[i] = INFINITY; }
  }
  //remove bad diagonals
  const std::vector<int> bad_diags = as<std::vector<int> >(outliers["bad.diagonals"]);
  const std::set<int> bad_diags_set(bad_diags.begin(), bad_diags.end());
  for (unsigned i=0; i<nbetas; ++i) {
    int diag_idx = bin2[i]-bin1[i];
    if (bad_diags_set.find(diag_idx) != bad_diags_set.end()) phihat_var[i] = INFINITY;
  }
}

void cts_to_signal_mat(const SignalRawData& raw, double eCprime, const Rcpp::NumericVector& beta_phi, SignalBinnedData& binned) {
    //extract data from holder
    const DataFrame& cts = raw.get_cts();
    int nbins = raw.get_nbins();
    double dispersion = raw.get_dispersion();
    const std::vector<double> phi = Rcpp::as<std::vector<double> >(beta_phi);
    const List& outliers = raw.get_outliers();

    //inputs
    int N = cts.nrows();
    std::vector<int> cts_bin1 = as<std::vector<int> >(cts["bin1"]);
    std::vector<int> cts_bin2 = as<std::vector<int> >(cts["bin2"]);
    std::vector<double> count = as<std::vector<double> >(cts["count"]);
    std::vector<double> lmu_nosig = as<std::vector<double> >(cts["lmu.nosig"]);
    std::vector<double> weight = as<std::vector<double> >(cts["weight"]); //not related to 1/var !
    
    //outputs
    int nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<double> phihat(nbetas, 0); //vectorized form
    std::vector<double> phihat_var(nbetas, 0);
    std::vector<double> ncounts(nbetas, 0);
    std::vector<int> bin1(nbetas, 0);
    std::vector<int> bin2(nbetas, 0);

    //build mat from cts
    double* pphi = const_cast<double*>(&phi[0]);
    cts_to_signal_mat_core(N, &cts_bin1[0], &cts_bin2[0], &count[0], &lmu_nosig[0],
                           &weight[0], nbins, dispersion,
                           pphi, eCprime, &phihat[0], &phihat_var[0], &ncounts[0], &bin1[0], &bin2[0]);
    //remove outliers
    remove_outliers(bin1, bin2, phihat_var, outliers);
                           
    //report back
    IntegerVector bin1_i, bin2_i, didx_i, dgrp_i;
    NumericVector phihat_i, phihat_var_i, ncounts_i, weight_i;
    bin1_i = wrap(bin1);
    bin2_i = wrap(bin2);
    bin1_i.attr("levels") = as<IntegerVector>(cts["bin1"]).attr("levels");
    bin2_i.attr("levels") = as<IntegerVector>(cts["bin2"]).attr("levels");
    bin1_i.attr("class") = CharacterVector::create("ordered", "factor");
    bin2_i.attr("class") = CharacterVector::create("ordered", "factor");
    phihat_i = wrap(phihat);
    phihat_var_i = wrap(phihat_var);
    ncounts_i = wrap(ncounts);
    weight_i = 1/phihat_var_i;
    didx_i = bin2_i-bin1_i;
    dgrp_i = floor(log10(bin2_i-bin1_i+1)*3);

    binned.set_bin1(bin1_i);
    binned.set_bin2(bin2_i);
    binned.set_phihat(phihat_i);
    binned.set_weight(weight_i);
    binned.set_ncounts(ncounts_i);
    binned.set_diag_idx(didx_i);
    binned.set_diag_grp(dgrp_i);
}

void cts_to_diff_mat(const DifferenceRawData& raw, const Rcpp::NumericVector& phi_ref, const Rcpp::NumericVector& beta_delta, DifferenceBinnedData& binned) {
    //assume eCprime = 0 for difference step
    const double eCprime=0;
    //compute ref matrix
    SignalRawData sraw_ref(raw.get_nbins(), raw.get_dispersion(), raw.get_ref(), raw.get_outliers());
    SignalBinnedData sbinned_ref;
    cts_to_signal_mat(sraw_ref, eCprime, phi_ref, sbinned_ref);
    //compute other matrix
    Rcpp::NumericVector phi_oth = phi_ref+beta_delta; //unthresholded
    SignalRawData sraw_oth(raw.get_nbins(), raw.get_dispersion(), raw.get_cts(), raw.get_outliers());
    SignalBinnedData sbinned_oth;
    cts_to_signal_mat(sraw_oth, eCprime, phi_oth, sbinned_oth);
    //report
    binned.set_bin1(sbinned_ref.get_bin1());
    binned.set_bin2(sbinned_ref.get_bin2());
    binned.set_beta_delta(beta_delta);
    binned.set_weight(sbinned_oth.get_weight());
    binned.set_deltahat(sbinned_oth.get_phihat() - sbinned_ref.get_phihat());
    binned.set_ncounts(sbinned_oth.get_ncounts() + sbinned_ref.get_ncounts());
    binned.set_diag_idx(sbinned_ref.get_diag_idx());
    binned.set_diag_grp(sbinned_ref.get_diag_grp());
    binned.set_weight_ref(sbinned_ref.get_weight());
    binned.set_phihat_ref(sbinned_ref.get_phihat());
    binned.set_phi_ref(phi_ref);
}

