#include <Rcpp.h>
using namespace Rcpp;

#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"

RCPP_MODULE(csnorm_cpp) {
    using namespace Rcpp ;

    function("wgfl_signal_BIC", &wgfl_signal_BIC, "documentation for wgfl_signal_BIC ");
    
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
}

