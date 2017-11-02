#include <Rcpp.h>
using namespace Rcpp;

#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "fast_binless.hpp"

RCPP_MODULE(binless_cpp) {
    using namespace Rcpp ;

    function("wgfl_signal_BIC", &wgfl_signal_BIC, "documentation for wgfl_signal_BIC ");
    
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
    
    function("fast_binless", &fast_binless, "documentation for fast_binless ");
    function("fast_binless_eval_cv", &fast_binless_eval_cv, "documentation for fast_binless_eval_cv ");
    function("fast_binless_difference", &fast_binless_difference, "documentation for fast_binless_difference ");
}

