#include <RcppEigen.h>
using namespace Rcpp;

#include "perf_iteration_signal.hpp"
#include "perf_iteration_diff.hpp"
#include "fast_binless.hpp"
#include "spline.hpp"
#include "cts_to_mat.hpp"

#include "FastData.hpp"
#include "fast_binless.hpp"

RCPP_EXPOSED_CLASS_NODECL(binless::fast::FastSignalData);
RCPP_EXPOSED_CLASS_NODECL(binless::fast::DecayConfig);


RCPP_MODULE(binless_cpp) {
    using namespace Rcpp ;

    function("wgfl_signal_BIC", &wgfl_signal_BIC, "documentation for wgfl_signal_BIC ");
    
    function("wgfl_diff_BIC", &wgfl_diff_BIC, "documentation for wgfl_diff_BIC ");
    
    function("fast_binless", &binless::fast::binless, "documentation for fast_binless ");
    function("fast_binless_eval_cv", &binless::fast::binless_eval_cv, "documentation for fast_binless_eval_cv ");
    function("fast_binless_difference", &binless::fast::binless_difference, "documentation for fast_binless_difference ");
    
    function("generate_spline_base", &generate_spline_base, "documentation for generate_spline_base ");
    
    function("rcpp_cts_to_signal_mat", &rcpp_cts_to_signal_mat, "documentation for rcpp_cts_to_signal_mat ");
    
    class_<binless::fast::DecayConfig>("DecayConfig")
     .constructor<double,double>()
     .field("Kdiag", &binless::fast::DecayConfig::Kdiag)
     .field("max_iter", &binless::fast::DecayConfig::max_iter)
     .field("sigma", &binless::fast::DecayConfig::sigma)
     .field("bins_per_bf", &binless::fast::DecayConfig::bins_per_bf)
     .field("tol_val", &binless::fast::DecayConfig::tol_val)
     .field("free_decay", &binless::fast::DecayConfig::free_decay)
    ;
    
    class_<binless::fast::DecayEstimator>("DecayEstimator")
     .constructor<const binless::fast::FastSignalData&, const binless::fast::DecayConfig&>()
    ;
}

