#include <Rcpp.h>
using namespace Rcpp;
#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "graph_trails.hpp"


RCPP_MODULE(gfl){
  using namespace Rcpp ;
  
  function("cts_to_signal_mat" , &cts_to_signal_mat  , "documentation for cts_to_signal_mat ");
  function("wgfl_signal_perf_warm" , &wgfl_signal_perf_warm  , "documentation for wgfl_signal_perf_warm ");
  
  function("cts_to_diff_mat" , &cts_to_diff_mat  , "documentation for cts_to_diff_mat ");
  function("wgfl_diff_perf_warm" , &wgfl_diff_perf_warm  , "documentation for wgfl_diff_perf_warm ");
  
  function("boost_triangle_grid_chain", &boost_triangle_grid_chain, "documentation for boost_triangle_grid_chain ");
  function("boost_chains_to_trails", &boost_chains_to_trails, "documentation for boost_chains_to_trails ");
  function("boost_build_patch_graph_components" , &boost_build_patch_graph_components  , "documentation for boost_build_patch_graph_components ");
  
} 

