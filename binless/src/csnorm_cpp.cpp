#include <RcppEigen.h>
using namespace Rcpp;

#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "fast_binless.hpp"
#include "spline.hpp"

RCPP_MODULE(binless_cpp) {
    using namespace Rcpp ;

    function("wgfl_signal_BIC", &wgfl_signal_BIC, "documentation for wgfl_signal_BIC ");
    
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
    
    function("fast_binless", &binless::fast::binless, "documentation for fast_binless ");
    function("fast_binless_eval_cv", &binless::fast::binless_eval_cv, "documentation for fast_binless_eval_cv ");
    function("fast_binless_difference", &binless::fast::binless_difference, "documentation for fast_binless_difference ");
    
    function("generate_spline_base", &generate_spline_base, "documentation for generate_spline_base ");
}

