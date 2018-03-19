#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <set>

#include "cts_to_mat.hpp"

#include "cts_core.h" //cts_to_signal_mat_core

#include "RawData.hpp"
#include "BinnedData.hpp"
#include "Traits.hpp"


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

void cts_to_signal_mat(const RawData<Signal>& raw, BinnedData<Signal>& binned) {
    //extract data from holder
    const DataFrame& cts = raw.get_cts();
    int nbins = raw.get_nbins();
    double dispersion = raw.get_dispersion();
    const std::vector<double> phi = Rcpp::as<std::vector<double> >(binned.get_beta_phi());
    const List& metadata = raw.get_metadata();

    //inputs
    int N = cts.nrows();
    std::vector<int> cts_bin1 = as<std::vector<int> >(cts["bin1"]);
    std::vector<int> cts_bin2 = as<std::vector<int> >(cts["bin2"]);
    std::vector<double> count = as<std::vector<double> >(cts["count"]);
    std::vector<double> lmu_nosig = as<std::vector<double> >(cts["lmu.nosig"]);
    std::vector<double> nobs = as<std::vector<double> >(cts["nobs"]);
    std::vector<double> log_decay = as<std::vector<double> >(cts["log_decay"]);
    
    //outputs
    int nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<double> phihat(nbetas, 0); //vectorized form
    std::vector<double> phihat_var(nbetas, 0);
    std::vector<double> phihat_var_nodecay(nbetas, 0);
    std::vector<double> ncounts(nbetas, 0);
    std::vector<int> bin1(nbetas, 0);
    std::vector<int> bin2(nbetas, 0);

    //build mat from cts
    double* pphi = const_cast<double*>(&phi[0]);
    cts_to_signal_mat_core(N, &cts_bin1[0], &cts_bin2[0], &count[0], &lmu_nosig[0],
                           &nobs[0], &log_decay[0], nbins, dispersion,
                           pphi, &phihat[0], &phihat_var[0], &phihat_var_nodecay[0], &ncounts[0], &bin1[0], &bin2[0]);
    //remove outliers
    remove_outliers(bin1, bin2, phihat_var, metadata);
                           
    //report back
    IntegerVector bin1_i, bin2_i, didx_i;
    NumericVector phihat_i, phihat_var_i, phihat_var_nodecay_i, ncounts_i, weight_i, weight_nodecay_i;
    bin1_i = wrap(bin1);
    bin2_i = wrap(bin2);
    bin1_i.attr("levels") = as<IntegerVector>(cts["bin1"]).attr("levels");
    bin2_i.attr("levels") = as<IntegerVector>(cts["bin2"]).attr("levels");
    bin1_i.attr("class") = CharacterVector::create("ordered", "factor");
    bin2_i.attr("class") = CharacterVector::create("ordered", "factor");
    phihat_i = wrap(phihat);
    phihat_var_i = wrap(phihat_var);
    phihat_var_nodecay_i = wrap(phihat_var_nodecay);
    ncounts_i = wrap(ncounts);
    weight_i = 1/phihat_var_i;
    weight_nodecay_i = 1/phihat_var_nodecay_i;
    didx_i = bin2_i-bin1_i;
    IntegerVector diag_grps = metadata["diag.grp"];
    IntegerVector dgrp_i(didx_i.size());
    for (unsigned i=0; i<didx_i.size(); ++i) {
        dgrp_i(i) = diag_grps(didx_i(i));
    }

    binned.set_bin1(bin1_i);
    binned.set_bin2(bin2_i);
    binned.set_phihat(phihat_i);
    binned.set_weight(weight_i);
    binned.set_weight_nodecay(weight_nodecay_i);
    binned.set_nobs(ncounts_i);
    binned.set_diag_idx(didx_i);
    binned.set_diag_grp(dgrp_i);
}

void cts_to_diff_mat(const RawData<Difference>& raw, BinnedData<Difference>& binned) {
    //compute ref matrix
    RawData<Signal> sraw_ref(raw.get_nbins(), raw.get_dispersion(), raw.get_ref(), raw.get_metadata());
    BinnedData<Signal> sbinned_ref;
    sbinned_ref.set_beta_phi(binned.get_phi_ref());
    cts_to_signal_mat(sraw_ref, sbinned_ref);
    //compute other matrix
    Rcpp::NumericVector phi_oth = binned.get_phi_ref()+binned.get_beta_delta(); //unthresholded
    RawData<Signal> sraw_oth(raw.get_nbins(), raw.get_dispersion(), raw.get_cts(), raw.get_metadata());
    BinnedData<Signal> sbinned_oth;
    sbinned_oth.set_beta_phi(phi_oth);
    cts_to_signal_mat(sraw_oth, sbinned_oth);
    //report
    binned.set_bin1(sbinned_ref.get_bin1());
    binned.set_bin2(sbinned_ref.get_bin2());
    binned.set_weight(sbinned_oth.get_weight());
    binned.set_weight_nodecay(sbinned_oth.get_weight_nodecay());
    binned.set_deltahat(sbinned_oth.get_phihat() - binned.get_phi_ref());
    binned.set_nobs(sbinned_oth.get_nobs() + sbinned_ref.get_nobs());
    binned.set_diag_idx(sbinned_ref.get_diag_idx());
    binned.set_diag_grp(sbinned_ref.get_diag_grp());
    binned.set_weight_ref(sbinned_ref.get_weight());
    binned.set_weight_nodecay_ref(sbinned_ref.get_weight_nodecay());
    binned.set_phihat_ref(sbinned_ref.get_phihat());
}

Rcpp::DataFrame rcpp_cts_to_signal_mat(int nbins, double alpha, const Rcpp::DataFrame cts, const Rcpp::List metadata) {
  RawData<Signal> raw(nbins, alpha, cts, metadata);
  BinnedData<Signal> binned;
  Rcpp::NumericVector beta(nbins*(nbins+1)/2);
  binned.set_beta_phi(beta);
  cts_to_signal_mat(raw,binned);
  binned.get_bin1();
  binned.get_nobs();
  binned.get_weight();
  binned.get_diag_idx();
  binned.get_diag_grp();
  binned.get_beta();
  DataFrame mat = DataFrame::create(_["bin1"]=binned.get_bin1(),
                                    _["bin2"]=binned.get_bin2(),
                                    _["phihat"]=binned.get_phihat(),
                                    _["nobs"]=binned.get_nobs(),
                                    _["weight"]=binned.get_weight(),
                                    _["diag.idx"]=binned.get_diag_idx(),
                                    _["diag.grp"]=binned.get_diag_grp(),
                                    _["beta"]=binned.get_beta());
  return mat;
}


