#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include "util.hpp"

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1) {
    std::vector<double> phi;
    phi.reserve(beta.size());
    for (std::vector<double>::const_iterator it = beta.begin(); it != beta.end();
            ++it) {
        double val = *it - eCprime;
        phi.push_back((val > 0) ? std::max(0., val - lam1) : std::min(0., val + lam1));
    }
    return phi;
}

Rcpp::NumericVector compute_phi_ref(const BinnedData<Difference>& binned,
                                    const Rcpp::NumericVector& delta) {
  const int N = delta.size();
  std::vector<double> phihat = Rcpp::as<std::vector<double> >(binned.get_phihat());
  std::vector<double> phihat_var = Rcpp::as<std::vector<double> >(binned.get_phihat_var());
  std::vector<double> phihat_ref = Rcpp::as<std::vector<double> >(binned.get_phihat_ref());
  std::vector<double> phihat_var_ref = Rcpp::as<std::vector<double> >(binned.get_phihat_var_ref());
  std::vector<double> phi_ref_r;
  phi_ref_r.reserve(N);
  for (int i=0; i<N; ++i) {
    double val;
    if (phihat_var_ref[i]==INFINITY && phihat_var[i]==INFINITY) {
      val=(phihat_ref[i]+phihat[i])/2;
    } else {
      val=(phihat_ref[i]/phihat_var_ref[i] + (phihat[i]-delta[i])/phihat_var[i])
      /(1/phihat_var_ref[i]+1/phihat_var[i]);
    }
    val = std::max(val,0.);
    phi_ref_r.push_back(val);
  }
  return Rcpp::wrap(phi_ref_r);
}

double get_maximum_admissible_weight(const Rcpp::NumericVector& w, double percentile) {
    //get sorted values of weights
    std::vector<double> sw = Rcpp::as<std::vector<double> >(w);
    std::sort(sw.begin(),sw.end());
    //find first nonzero weight
    unsigned i;
    for (i=0; i<sw.size(); ++i) if (sw[i]>0) break;
    //find percentile of positive values
    i = unsigned(percentile*(sw.size()-1-i)/100.)+i;
    //return its value
    return sw[i];
}

BinDesign make_even_bins(double lower, double upper, unsigned Nbins) {
  double size = upper - lower;
  std::map<double,unsigned> ret;
  std::vector<double> begins, ends;
  for (unsigned i=0; i<=Nbins; ++i) { //size of map is Nbins+1
    double begin = lower + i*size/double(Nbins);
    ret[begin]=i;
    if (i<Nbins) begins.push_back(begin);
    if (i>0) ends.push_back(begin);
  }
  return BinDesign{Nbins,ret,begins,ends};
}

Eigen::SparseMatrix<double, Eigen::RowMajor> bin_data(const Eigen::VectorXd& data, const BinDesign& design, bool drop) {
  unsigned Ndata = data.rows();
  std::vector<Eigen::Triplet<double> > tripletList;
  tripletList.reserve(Ndata);
  for (unsigned i=0; i<Ndata; ++i) {
    auto it = design.map.upper_bound(data(i));
    unsigned j;
    if (it == design.map.end()) {
      j = design.map.rbegin()->second -1; //place in last bin
    } else {
      j = it->second - 1; //upper_bound returns pointer to next bin
    }
    tripletList.push_back(Eigen::Triplet<double>(j,i,1.));
  }
  Eigen::SparseMatrix<double, Eigen::RowMajor> mat(design.Nbins,Ndata);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  if (drop) {
    Eigen::SparseMatrix<double, Eigen::RowMajor> ret(design.Nbins,Ndata);
    unsigned Nrow=0;
    for (unsigned i=0; i<design.Nbins; ++i) {
      auto mrow = mat.row(i);
      if (mrow.sum()>0) {
        ret.row(Nrow++) = mrow;
      }
    }
    ret.conservativeResize(Nrow,Ndata);
    ret.makeCompressed();
    return ret;
  } else {
    mat.makeCompressed();
    return mat;
  }
}

Eigen::SparseMatrix<double, Eigen::RowMajor> bin_data_evenly(const Eigen::VectorXd& data, unsigned Nbins, bool drop) {
  auto design = make_even_bins(data.minCoeff(),data.maxCoeff(),Nbins);
  return bin_data(data,design,drop);
}

Rcpp::DataFrame create_empty_matrix_cpp(Rcpp::IntegerVector name_ordered, Rcpp::IntegerVector bins_ordered) {
  //extract info and create vectors
  Rcpp::CharacterVector name_labels(name_ordered.attr("levels"));
  Rcpp::CharacterVector bin_labels(bins_ordered.attr("levels"));
  unsigned n_names = name_labels.size();
  unsigned n_bins = bin_labels.size();
  std::vector<int> name,bin1,bin2;
  unsigned nrows = n_names * (n_bins * (n_bins+1))/2;
  name.reserve(nrows);
  bin1.reserve(nrows);
  bin2.reserve(nrows);
  //build vectors
  for (int n=1; n <= n_names; ++n) {
    for (int i=1; i <= n_bins; ++i) {
      for (int j=i; j <= n_bins; ++j) {
        name.push_back(n);
        bin1.push_back(i);
        bin2.push_back(j);
      }
    }
  }
  //convert to rcpp
  Rcpp::IntegerVector rname(Rcpp::wrap(name));
  rname.attr("levels") = name_labels;
  rname.attr("class") = CharacterVector::create("ordered", "factor");
  Rcpp::IntegerVector rbin1(Rcpp::wrap(bin1));
  rbin1.attr("levels") = bin_labels;
  rbin1.attr("class") = CharacterVector::create("ordered", "factor");
  Rcpp::IntegerVector rbin2(Rcpp::wrap(bin2));
  rbin2.attr("levels") = bin_labels;
  rbin2.attr("class") = CharacterVector::create("ordered", "factor");
  //return data frame
  return Rcpp::DataFrame::create(_["name"]=rname, _["bin1"]=rbin1, _["bin2"]=rbin2);
}





