#include <Rcpp.h>
using namespace Rcpp;

#include "cts_to_mat.hpp"
#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "graph_helpers.hpp"

RCPP_MODULE(csnorm_cpp) {
    using namespace Rcpp ;

    function("wgfl_signal_perf_warm", &wgfl_signal_perf_warm,
             "documentation for wgfl_signal_perf_warm ");
    function("wgfl_signal_BIC", &wgfl_signal_BIC,
             "documentation for wgfl_signal_BIC ");
    function("wgfl_signal_BIC_fixed", &wgfl_signal_BIC_fixed,
             "documentation for wgfl_signal_BIC_fixed ");

    function("wgfl_diff_perf_warm", &wgfl_diff_perf_warm,
             "documentation for wgfl_diff_perf_warm ");
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
    function("wgfl_diff_BIC_fixed", &wgfl_diff_BIC_fixed, "documentation for wgfl_diff_BIC_fixed ");
    
    function("get_patch_numbers",
             &get_patch_numbers,
             "documentation for get_patch_numbers ");
}

