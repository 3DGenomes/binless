#include <Rcpp.h>
#include <vector>

#include "SignalWeightsUpdater.hpp"
#include "cts_to_mat.hpp"

void SignalWeightsUpdater::update(const std::vector<double>& beta) {
    mat_ = cts_to_signal_mat(cts_, nrows_, dispersion_, beta, 0, outliers_);
}
