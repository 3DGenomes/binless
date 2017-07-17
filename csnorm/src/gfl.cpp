#include <Rcpp.h>
using namespace Rcpp;

#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "graph_trails.hpp"

RCPP_MODULE(gfl) {
    using namespace Rcpp ;

    function("cts_to_signal_mat", &cts_to_signal_mat,
             "documentation for cts_to_signal_mat ");
    function("wgfl_signal_perf_warm", &wgfl_signal_perf_warm,
             "documentation for wgfl_signal_perf_warm ");
    function("wgfl_signal_perf_opt_lambda1_eCprime",
             &wgfl_signal_perf_opt_lambda1_eCprime,
             "documentation for wgfl_signal_perf_opt_lambda1_eCprime");
    function("wgfl_signal_BIC", &wgfl_signal_BIC,
             "documentation for wgfl_signal_BIC ");
    function("wgfl_signal_BIC_fixed", &wgfl_signal_BIC_fixed,
             "documentation for wgfl_signal_BIC_fixed ");

    function("cts_to_diff_mat", &cts_to_diff_mat,
             "documentation for cts_to_diff_mat ");
    function("wgfl_diff_perf_warm", &wgfl_diff_perf_warm,
             "documentation for wgfl_diff_perf_warm ");
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
    function("wgfl_diff_BIC_fixed", &wgfl_diff_BIC_fixed, "documentation for wgfl_diff_BIC_fixed ");
    
    function("boost_triangle_grid_chain", &boost_triangle_grid_chain,
             "documentation for boost_triangle_grid_chain ");
    function("boost_chains_to_trails", &boost_chains_to_trails,
             "documentation for boost_chains_to_trails ");
    function("boost_build_patch_graph_components",
             &boost_build_patch_graph_components,
             "documentation for boost_build_patch_graph_components ");
}

