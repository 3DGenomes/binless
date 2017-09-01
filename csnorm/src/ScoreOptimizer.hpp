#ifndef SCORE_OPTIMIZER_HPP
#define SCORE_OPTIMIZER_HPP

#include <Rcpp.h>
#include <vector>

//Find optimum value based on a list of score evaluations
template<typename Score> class ScoreOptimizer {
public:
    typedef Rcpp::NumericVector value_t;
    
    value_t optimize(const std::vector<value_t>& values) const {
        if (values.size()==0) Rcpp::stop("Found 0 values!\n");
        //find minimum
        auto minval_it = std::min_element(values.cbegin(), values.cend(),
                    [](const value_t& a, const value_t& b) {return Rcpp::is_true(a["BIC"]<b["BIC"]); } );
        //optionally, go upwards
        return best = add_kSD<Score::k>(minval_it, values.cend());
    }
    
private:
    template<unsigned k> value_t add_kSD(const std::vector<value_t>::const_iterator minval_it,
                                         const std::vector<value_t>::const_iterator minval_end) const {
        auto it = minval_it;
        const double minscore = (*it)["BIC"] + k * (*it)["BIC.sd"];
        do { ++it; } while (it != minval_end && Rcpp::as<double>((*it)["BIC"]) <= minscore );
        if (it == minval_end) --it;
        return *it;
    }
    
    template<> value_t add_kSD<0>(std::vector<value_t>::const_iterator minval_it,
                                  const std::vector<value_t>::const_iterator) const {
        return (*minval_it);
    }
    
};
#endif
