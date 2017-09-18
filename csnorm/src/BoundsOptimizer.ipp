#include <Rcpp.h>

#include "Traits.hpp"

template<typename Sign>
bounds_t BoundsOptimizer<EstimatedOffset, Sign>::optimize_bounds(double val) const {
    //split data in two groups
    LogicalVector grp1 = beta_ <= val;
    Rcpp::NumericVector w1 = weight_[grp1];
    Rcpp::NumericVector w2 = weight_[!grp1];
    bool grp2empty = w2.size() == 0 || sum(w2) == 0;
    bool grp1empty = w1.size() == 0 || sum(w1) == 0;
    double UB,LB;
    
    //compute UB and LB
    if ((!grp2empty) && (!grp1empty)) {//non-degenerate case
        Rcpp::NumericVector x1 = beta_[grp1];
        UB = max(x1);
        LB = xmin_;
    } else if (grp2empty) {//grp2 is empty, UB > xmax_
        UB = xmax_;
        LB = xmin_;
    } else {//grp1 is empty, UB <= xmin_
        UB = xmin_;
        LB = xmin_;
    }
    if (UB<LB) Rcpp::stop("Found LB > UB!\n");
    return bounds_t{LB, UB};
}


template<typename Sign>
bounds_t BoundsOptimizer<ZeroOffset, Sign>::optimize_bounds(double val) const {
    //split data in two groups and determine constraint betas
    Rcpp::LogicalVector grp1 = beta_ > val, grp2 = beta_ < -val;
    Rcpp::NumericVector w1 = weight_[grp1];
    Rcpp::NumericVector w2 = weight_[grp2];
    bool grp1empty = w1.size() == 0 || sum(w1) == 0;
    bool grp2empty = w2.size() == 0 || sum(w2) == 0;
    double xk,xkp1;
    if (val>maxabsval_) {
        xk=maxabsval_;
        xkp1=std::numeric_limits<double>::infinity();
    } else if (val<minabsval_) {
        xk=0;
        xkp1=minabsval_;
    } else {
        xk=max(Rcpp::as<Rcpp::NumericVector>(absval_[absval_<=val]));
        xkp1=min(Rcpp::as<Rcpp::NumericVector>(absval_[absval_>val]));
    }
    double UB = xk;
    return bounds_t{-UB,UB};
}

