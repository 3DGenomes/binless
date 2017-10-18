#include <Rcpp.h>
using namespace Rcpp;

#include "fast_binless.hpp"

RCPP_MODULE(binless_cpp) {
    using namespace Rcpp ;

    function("fast_binless", &fast_binless, "documentation for fast_binless ");
    function("fast_binless_eval_cv", &fast_binless_eval_cv, "documentation for fast_binless_eval_cv ");
    function("fast_binless_difference", &fast_binless_difference, "documentation for fast_binless_difference ");
}

