#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <set>

#include "cts_to_mat.hpp"

#include "cts_core.h" //cts_to_signal_mat_core

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

DataFrame cts_to_signal_mat(const DataFrame& cts, int nbins, double dispersion,
                            const std::vector<double>& phi,
                            double eCprime, List outliers) {
    //inputs
    int N = cts.nrows();
    std::vector<int> cts_bin1 = as<std::vector<int> >(cts["bin1"]);
    std::vector<int> cts_bin2 = as<std::vector<int> >(cts["bin2"]);
    std::vector<double> count = as<std::vector<double> >(cts["count"]);
    std::vector<double> lmu_nosig = as<std::vector<double> >(cts["lmu.nosig"]);
    std::vector<double> weight = as<std::vector<double> >(cts["weight"]);
    
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

    
    return DataFrame::create(_["bin1"]=bin1_i, _["bin2"]=bin2_i,
                             _["phihat"]=phihat_i,
                             _["phihat.var"]=phihat_var_i, _["ncounts"]=ncounts_i, _["weight"]=weight_i,
                             _["diag.idx"]=didx_i, _["diag.grp"]=dgrp_i);
}

DataFrame cts_to_diff_mat(const DataFrame& cts, const DataFrame ref, int nbins,
                          double dispersion,
                          const std::vector<double>& phi_ref, const std::vector<double>& delta, List outliers) {
    //assume eCprime = 0 for difference step
    const double eCprime=0;
    const DataFrame mat_ref = cts_to_signal_mat(ref, nbins, dispersion, phi_ref,
                                                eCprime, outliers);
    
    std::vector<double> phi_oth;
    phi_oth.reserve(delta.size());
    for (int i=0; i<delta.size(); ++i) phi_oth[i] = phi_ref[i] + delta[i];
    const DataFrame mat_oth = cts_to_signal_mat(cts, nbins, dispersion, phi_oth,
                                                eCprime, outliers);
    IntegerVector bin1 = mat_ref["bin1"];
    IntegerVector bin2 = mat_ref["bin2"];
    NumericVector phihat_ref = mat_ref["phihat"];
    NumericVector phihat_var_ref = mat_ref["phihat.var"];
    NumericVector phihat = mat_oth["phihat"];
    NumericVector phihat_var = mat_oth["phihat.var"];
    NumericVector ncounts_oth = mat_oth["ncounts"];
    NumericVector ncounts_ref = mat_ref["ncounts"];
    NumericVector deltahat = phihat-phihat_ref;
    NumericVector deltahat_var = phihat_var+phihat_var_ref;
    NumericVector ncounts = ncounts_ref+ncounts_oth;
    NumericVector weight = 1/deltahat_var;
    IntegerVector didx = mat_ref["diag.idx"];
    IntegerVector dgrp = mat_ref["diag.grp"];
    return DataFrame::create(_["bin1"]=bin1, _["bin2"]=bin2,
                             _["phihat"]=phihat, _["phihat.var"]=phihat_var,
                             _["phihat.ref"]=phihat_ref, _["phihat.var.ref"]=phihat_var_ref,
                             _["deltahat"]=deltahat, _["deltahat.var"]=deltahat_var, _["ncounts"]=ncounts,
                             _["weight"]=weight, _["diag.idx"]=didx, _["diag.grp"]=dgrp);
}

