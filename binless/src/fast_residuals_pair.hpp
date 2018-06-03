#ifndef FAST_RESIDUALS_PAIR_HPP
#define FAST_RESIDUALS_PAIR_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

namespace binless {
namespace fast {

struct ResidualsPair { std::vector<double> residuals,weights; };

}
}

#endif

