#include <Rcpp.h>

#include "Tags.hpp"

template<typename Sign>
bounds_t BoundsComputer<EstimatedOffset, Sign>::optimize_bounds(double val) const {
    //split data in two groups
    LogicalVector grp1 = value_ <= val;
    Rcpp::NumericVector w1 = weight_[grp1];
    Rcpp::NumericVector w2 = weight_[!grp1];
    bool grp2empty = w2.size() == 0 || sum(w2) == 0;
    bool grp1empty = w1.size() == 0 || sum(w1) == 0;
    double b1,b2, xk, xkp1;
    double a,b,UB,LB;
    double xmin = min(value_), xmax = max(value_); //here, we assume xmin >= minval_
    
    //compute UB and LB
    if ((!grp2empty) && (!grp1empty)) {//non-degenerate case
        Rcpp::NumericVector betahat1 = valuehat_[grp1];
        Rcpp::NumericVector betahat2 = valuehat_[!grp1];
        Rcpp::NumericVector x1 = value_[grp1];
        Rcpp::NumericVector x2 = value_[!grp1];
        b1 = sum(w1*betahat1)/sum(w1);
        b2 = sum(w2*(x2-betahat2))/sum(w2);
        xk = max(x1);
        xkp1 = min(x2);
        //compute unconstrained minimum for UB and LB
        a = b1 + b2;
        b = b1 - b2;
        //compute optimal UB and LB
        UB = std::max(std::min(xkp1, a), xk);
        LB = std::min(b, minval_);
    } else if (grp2empty) {//grp2 is empty, UB > xmax
        Rcpp::NumericVector w1 = weight_;
        Rcpp::NumericVector betahat1 = valuehat_;
        b1 = sum(w1*betahat1)/sum(w1);
        b2 = 0;
        xk = xmax;
        xkp1 = std::numeric_limits<double>::infinity();
        //compute unconstrained minimum for UB and LB
        a = b1 + b2;
        b = b1 - b2;
        //compute optimal UB and LB
        if (b1 >= (minval_+xmax)/2.) {
            LB = minval_; //or any value < minval_
            UB = 2*b1 - LB;
        } else {
            UB = xmax; //or any value > xmax
            LB = 2*b1 - UB;
        }
    } else {//grp1 is empty, UB <= xmin
        if (!grp1empty) Rcpp::Rcout << "This should not have happened!\n";
        Rcpp::NumericVector w2 = weight_;
        Rcpp::NumericVector betahat2 = valuehat_;
        Rcpp::NumericVector x2 = value_;
        b1 = 0;
        b2 = sum(w2*(x2-betahat2))/sum(w2);
        xk = -std::numeric_limits<double>::infinity();
        xkp1 = xmin;
        //compute unconstrained minimum for UB and LB
        a = b1 + b2;
        b = b1 - b2;
        //compute optimal UB and LB
        if (b2 >= (xmin-minval_)/2.) {
            UB = xmin;
            LB = UB - 2*b2;
        } else {
            LB = minval_;
            UB = LB + 2*b2;
        }
    }
    if (UB<LB) Rcpp::stop("Found LB > UB!\n");
    return bounds_t{LB, UB};
}


template<typename Sign>
bounds_t BoundsComputer<ZeroOffset, Sign>::optimize_bounds(double val) const {
    //split data in two groups and determine constraint values
    //Rcout << "  GET at " << val << " maxabsval_= " << maxabsval_ << std::endl;
    Rcpp::LogicalVector grp1 = value_ > val, grp2 = value_ < -val;
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
    //Rcout << "  xk= " << xk << " xkp1= " << xkp1 << std::endl;
    //determine unconstrained UB
    double UB;
    //Rcout << "  grp1empty= " << grp1empty << " grp2empty= " << grp2empty << std::endl;
    if (grp1empty && grp2empty) {
        //UB larger than any value
        UB=maxabsval_;
    } else if (grp1empty) {
        //UB larger than any positive value
        Rcpp::NumericVector betahat2 = valuehat_[grp2];
        Rcpp::NumericVector x2 = value_[grp2];
        UB = -sum(w2*(x2-betahat2))/sum(w2);
    } else if (grp2empty) {
        Rcpp::NumericVector betahat1 = valuehat_[grp1];
        Rcpp::NumericVector x1 = value_[grp1];
        UB = sum(w1*(x1-betahat1))/sum(w1);
    } else {
        Rcpp::NumericVector betahat1 = valuehat_[grp1];
        Rcpp::NumericVector x1 = value_[grp1];
        double sw1 = sum(w1);
        double b1 = sum(w1*(x1-betahat1))/sw1;
        Rcpp::NumericVector betahat2 = valuehat_[grp2];
        Rcpp::NumericVector x2 = value_[grp2];
        double sw2 = sum(w2);
        double b2 = sum(w2*(x2-betahat2))/sw2;
        UB = (sw1*b1-sw2*b2)/(sw1+sw2);
    }
    //apply constraint
    UB = std::min(std::max(UB,xk),xkp1);
    if (minUB_ <= xkp1 && minUB_ > xk) UB=std::max(minUB_,UB);
    //Rcout << "  UB= " << UB << " minUB= " << minUB_ << std::endl;
    return bounds_t{-UB,UB};
}

