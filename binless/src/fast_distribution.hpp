#ifndef FAST_DISTRIBUTION_HPP
#define FAST_DISTRIBUTION_HPP

#include <Rcpp.h>
using namespace Rcpp;

namespace binless {
namespace fast {

//distribution families
struct Distribution {};

struct NormalDistribution : public Distribution {
  double sigma;
};

struct PoissonDistribution : public Distribution {};

struct NegativeBinomialDistribution : public Distribution {
  double alpha; // dispersion parameter (variance is mu + mu^2/alpha)
};


//methods to sample remaining parameters of distributions
// any Sampler specialization must have three calls
// - a constructor that stores a reference to the distribution
// - one init call to set the parameter to a sensible value
// - one sample call, to update their value
template<typename Dist> class Sampler {};

template<> class Sampler<PoissonDistribution> {
public:
  Sampler(PoissonDistribution&) {}
  void init() {}
  void sample() {}
};

template<> class Sampler<NegativeBinomialDistribution> {
public:
  Sampler(NegativeBinomialDistribution& nb) : nb_(nb) {}
  void init() { nb_.alpha = 1; }
  void sample() {}
private:
  NegativeBinomialDistribution& nb_;
};

}
}

#endif

