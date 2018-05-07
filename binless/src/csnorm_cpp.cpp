#include <RcppEigen.h>
using namespace Rcpp;

#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "fast_binless.hpp"
#include "spline.hpp"
#include "cts_to_mat.hpp"

RCPP_MODULE(binless_cpp) {
    using namespace Rcpp ;

    function("wgfl_signal_BIC", &wgfl_signal_BIC, "documentation for wgfl_signal_BIC ");
    
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
    
    function("fast_binless", &binless::fast::binless,
             List::create(_["obs"], _["nbins"], _["lam2"], _["alpha"], _["nouter"]=25, _["tol_val"]=2e-1, _["bg_steps"]=5, _["free_decay"]=10000),
             "documentation for fast_binless ");
    //function("fast_binless_eval_cv", &binless::fast::binless_eval_cv, "documentation for fast_binless_eval_cv ");
    function("fast_binless_difference", &binless::fast::binless_difference,
             List::create(_["obs"], _["lam2"], _["ref"], _["alpha"], _["tol_val"]=2e-1),
             "documentation for fast_binless_difference ");
    
    function("generate_spline_base", &generate_spline_base, "documentation for generate_spline_base ");
    
    function("rcpp_cts_to_signal_mat", &rcpp_cts_to_signal_mat, "documentation for rcpp_cts_to_signal_mat ");
}

