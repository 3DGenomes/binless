#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <string>

#include "optimize_lambda1_eCprime.hpp"
#include "util.hpp"
#include "graph_trails.hpp" //boost_build_patch_graph_components
#include <boost/math/tools/minima.hpp> //brent_find_minima

obj_lambda1_eCprime::obj_lambda1_eCprime(double minval, double maxval,
        double tol_val,
        bool constrained, IntegerVector patchno, NumericVector forbidden_vals,
        NumericVector value, NumericVector weight, NumericVector valuehat,
        NumericVector ncounts) : minval_(minval), valrange_(maxval-minval),
    tol_val_(tol_val), lsnc_(log(sum(ncounts))), constrained_(constrained),
    patchno_(patchno), forbidden_vals_(forbidden_vals),
    value_(value), weight_(weight), valuehat_(valuehat) {}

double obj_lambda1_eCprime::operator()(double val) const {
    return get(val, "opt")["BIC"];
    //return get(val)["BIC"];
}

//take val as a starting point for optimization of UB and LB at constant dof
NumericVector obj_lambda1_eCprime::get(double val, std::string msg) const {
    //split data in two groups
    LogicalVector grp1 = value_ <= val;
    double b1,b2, xk, xkp1;
    if (is_true(all(!grp1))) {
        //grp1 is empty, UB <= xmin
        b1 = 0;
        xk = -std::numeric_limits<double>::infinity();
        xkp1 = min(value_);
    } else {
        NumericVector w1 = weight_[grp1];
        NumericVector betahat1 = valuehat_[grp1];
        b1 = sum(w1*betahat1)/sum(w1);
    }
    if (is_true(all(grp1))) {
        //grp2 is empty, UB > xmax
        b2 = 0;
        xk = max(value_);
        xkp1 = std::numeric_limits<double>::infinity();
    } else {
        NumericVector w2 = weight_[!grp1];
        NumericVector betahat2 = valuehat_[!grp1];
        NumericVector x2 = value_[!grp1];
        b2 = sum(w2*(x2-betahat2))/sum(w2);
    }
    NumericVector x1 = value_[grp1];
    if (is_false(all(grp1)) && is_false(all(!grp1))) {
      NumericVector x2 = value_[!grp1];
      xk = max(x1);
      xkp1 = min(x2);
    }
    double xmin = min(value_);
    //compute unconstrained minimum for UB and LB
    double a = b1 + b2;
    double b = b1 - b2;
    //compute optimal UB and LB
    double UB = std::max(std::min(xkp1, a), xk);
    double LB = std::min(b, xmin);
    //compute eCprime, lambda1
    double eCprime = (UB+LB)/2;
    double lambda1 = (UB-LB)/2;
    //Rcout << "EVAL at val= " << val << " LB= " << LB << " UB= " << UB << " b1= " << b1 << " b2= " << b2 << " xmin= " << xmin << " xk= " << xk << " xkp1= " << xkp1 << std::endl; 
    //check if solution is feasible
    if (constrained_) {
        if ( is_true(any( (forbidden_vals_>UB+tol_val_/2) | (forbidden_vals_<LB-tol_val_/2) )) ) {
            if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1 << " eCprime= " << eCprime
                                    << " BIC= Inf dof= NA" << std::endl;
            return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                         _["BIC"]=std::numeric_limits<double>::max(), _["dof"]=NumericVector::get_na(),
                                         _["UB"]=UB, _["LB"]=LB);
        }
    }
    //compute dof and BIC
    std::vector<double> value_r = as<std::vector<double> >(value_);
    NumericVector soft = wrap(soft_threshold(value_r, eCprime, lambda1));
    IntegerVector selected = patchno_[abs(soft)>tol_val_/2];
    const int dof = unique(selected).size();
    const double BIC = sum(weight_ * SQUARE(valuehat_ - (soft + eCprime))) +
                       lsnc_*dof;
    if (!msg.empty()) Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= " << eCprime
                            << " BIC= " << BIC  << " dof= " << dof 
                            << " UB= " << UB  << " LB= " << LB << std::endl;
    return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                 _["BIC"]=BIC, _["dof"]=dof, _["UB"]=UB, _["LB"]=LB);
}

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num) {
    //extract vectors
    double lmin = std::max(lambda1_min,tol_val/2);
    std::clock_t c_start = std::clock();
    NumericVector weight = mat["weight"];
    NumericVector phihat = mat["phihat"];
    NumericVector beta = mat["beta"];
    NumericVector ncounts = mat["ncounts"];
    IntegerVector diag_idx = mat["diag.idx"];
    //get patch nos and sorted values
    List cl = boost_build_patch_graph_components(nbins, mat, tol_val);
    IntegerVector patchno = cl["membership"];
    NumericVector patchvals = get_patch_values(beta, patchno);
    //treat border case
    if (as<double>(cl["no"]) == 1) {
        Rcout << " OBJ final ok lambda1= " << lmin << " eCprime= " << patchvals(0)
                << " BIC= " << sum(weight * SQUARE(phihat)) << " dof= 0" << std::endl;
        return NumericVector::create(_["eCprime"]=min(patchvals)+tol_val/2, _["lambda1"]=lmin,
                                     _["dof"]=0,
                                     _["BIC"]=sum(weight * SQUARE(phihat)),
                                     _["c_init"]=-1, _["c_brent"]=-1, _["c_refine"]=-1);
    }
    double minval = patchvals(0);
    double maxval = patchvals(patchvals.size()-1);
    //if constraint is on, decay and signal must adjust so that
    //there is at least one zero signal value per diagonal idx
    NumericVector forbidden_vals;
    if (constrained) {
        forbidden_vals = get_minimum_diagonal_values(beta, diag_idx);
        lmin = std::max(lmin, (max(forbidden_vals)-minval)/2);
    }
    //create functor
    obj_lambda1_eCprime obj(minval, maxval, tol_val, constrained, patchno,
                            forbidden_vals,
                            beta, weight, phihat, ncounts);
    //for (int i=0; i<forbidden_vals.size(); ++i) Rcout << "fv[ " << i << " ]= "<< forbidden_vals[i] << std::endl;
    //loop over patch values
    std::clock_t c_in1 = std::clock();
    List val,best;
    best = obj.get(patchvals(0) - 2*tol_val, "opt"); //ensures  query is < xmin
    for (int i=0; i<patchvals.size(); ++i) {
      if (patchvals(i) <= max(forbidden_vals)) continue;
      val = obj.get(patchvals(i) + 2*tol_val, "opt");
      if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
    }
    val = obj.get(patchvals(patchvals.size()-1) + 2*tol_val, "opt"); //ensures  query is > xmax
    if (as<double>(val["BIC"]) < as<double>(best["BIC"])) best=val;
    std::clock_t c_in2 = std::clock();
    //finalize
    obj.get(as<double>(best["UB"])+2*tol_val,"final");
    return NumericVector::create(_["eCprime"]=best["eCprime"], _["lambda1"]=best["lambda1"],
                                 _["UB"]=best["UB"], _["LB"]=best["LB"],
                                 _["BIC"]=best["BIC"], _["dof"]=best["dof"],
                                 _["c_init"]=c_in1-c_start, _["c_brent"]=0, _["c_refine"]=c_in2-c_in1);
}
