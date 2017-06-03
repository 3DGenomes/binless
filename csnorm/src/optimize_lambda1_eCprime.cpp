#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>

#include "optimize_lambda1_eCprime.hpp"
#include "util.hpp" //soft_threshold and SQUARE
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include <boost/math/tools/minima.hpp> //brent_find_minima

obj_lambda1_eCprime::obj_lambda1_eCprime(double minval, double maxval, double tol_val,
                      bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
                      NumericVector value, NumericVector weight, NumericVector valuehat,
                      NumericVector ncounts) : minval_(minval), maxval_(maxval), valrange_(maxval-minval),
                      tol_val_(tol_val), lsnc_(log(sum(ncounts))), constrained_(constrained),
                      patchno_(patchno), forbidden_vals_(forbidden_vals),
                      value_(value), weight_(weight), valuehat_(valuehat) {}
  
double obj_lambda1_eCprime::operator()(double const& x) const { return get(std::pow(10,x))["BIC"]; }
  
NumericVector obj_lambda1_eCprime::get(double const& lambda1) const {
    double eCprime = (valrange_ < (2*lambda1+tol_val_)) ? (maxval_+minval_)/2. : lambda1+minval_-tol_val_;
    if (constrained_) {
      if (is_true(any(abs(eCprime-forbidden_vals_)>lambda1+tol_val_)))
        return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                     _["BIC"]=std::numeric_limits<double>::max(), _["dof"]=NumericVector::get_na());
    }
    std::vector<double> value_r = as<std::vector<double> >(value_);
    NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
    IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
    const int dof = unique(selected).size();
    const double BIC = sum(weight_ * SQUARE(valuehat_ - (soft + eCprime))) + lsnc_*dof;
    return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                 _["BIC"]=BIC, _["dof"]=dof);
}

NumericVector get_patch_values(NumericVector value, IntegerVector patchno) {
  //store value for each patch
  int npatches = max(patchno)+1;
  NumericVector unique_values(npatches); //patchno starts at 0
  for (int i=0; i<patchno.size(); ++i) unique_values(patchno(i)) = value(i);
  std::sort(unique_values.begin(), unique_values.end());
  return unique_values;
}

NumericVector get_minimum_diagonal_values(NumericVector value, IntegerVector diag_idx) {
  //store value for each patch
  int ndiags = max(diag_idx)+1;
  NumericVector diagvals(ndiags, max(value)); //diag_idx starts at 0
  for (int i=0; i<diag_idx.size(); ++i) diagvals(diag_idx(i)) = std::min(diagvals(diag_idx(i)),value(i));
  return diagvals;
}

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins, double tol_val, bool constrained, double lambda1_min) {
  //extract vectors
  NumericVector weight = mat["weight"];
  NumericVector phihat = mat["valuehat"];
  NumericVector phi = mat["value"];
  NumericVector ncounts = mat["ncounts"];
  IntegerVector diag_idx = mat["diag.idx"];
  //get patch nos and sorted values
  List cl = boost_build_patch_graph_components(nbins, mat, tol_val);
  IntegerVector patchno = cl["membership"];
  NumericVector patchvals = get_patch_values(phi, patchno);
  //treat border case
  if (as<double>(cl["no"]) == 1)
    return NumericVector::create(_["eCprime"]=patchvals(0), _["lambda1"]=lambda1_min, _["dof"]=0,
                        _["BIC"]=sum(weight * SQUARE(phihat - (phi + patchvals(0)))) );
  double minval = patchvals(0);
  double maxval = patchvals(patchvals.size()-1);
  //if constraint is on, decay and signal must adjust so that 
  //there is at least one zero signal value per diagonal idx
  NumericVector forbidden_vals;
  if (constrained) {
    forbidden_vals = get_minimum_diagonal_values(phi, diag_idx);
    lambda1_min = std::max(lambda1_min, (max(forbidden_vals)-minval)/2 + tol_val);
  }
  //create functor
  obj_lambda1_eCprime obj(minval, maxval, tol_val, constrained, patchno, forbidden_vals,
                          phi, weight, phihat, ncounts);
  //treat second border case
  if (maxval-minval <= 2*lambda1_min) return obj.get(lambda1_min);
  //optimize
  int bits = -8*std::log10(tol_val)+1;
  boost::uintmax_t maxiter = 100000;
  std::pair<double,double> ret = boost::math::tools::brent_find_minima(obj,
                                                   std::log10(std::max(lambda1_min,tol_val/2)),
                                                   std::log10(maxval-minval), bits, maxiter);
  return obj.get(pow(10,ret.first));
  //dt=foreach(lambda1=seq(lambda1.min,valrange,length.out=50), .combine=rbind) %do% obj(lambda1)
  //ggplot(dt)+geom_point(aes(lambda1,BIC))+geom_line(aes(lambda1,BIC))
}
