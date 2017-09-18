#include <Rcpp.h>

#include "Traits.hpp"

template<typename Sign>
bounds_t BoundsOptimizer<EstimatedOffset, Sign>::optimize_bounds(double val) const {
    Rcpp::Rcout << "\nBoundsOptimizer<EstimatedOffset, Sign> will optimize value " << val << "\n";
    //split data in two groups
    LogicalVector grp1 = beta_ <= val;
    Rcpp::NumericVector w1 = weight_[grp1];
    Rcpp::NumericVector w2 = weight_[!grp1];
    bool grp2empty = w2.size() == 0 || sum(w2) == 0;
    bool grp1empty = w1.size() == 0 || sum(w1) == 0;
    double b1,b2, xk, xkp1;
    double a,b,UB,LB;
    double xmin = min(beta_), xmax = max(beta_); //here, we assume xmin >= minval_
    Rcpp::Rcout << "grp1empty= " << grp1empty << " grp2empty= " << grp2empty << " xmin= " << xmin << " xmax= " << xmax << "\n";
    
    //compute UB and LB
    if ((!grp2empty) && (!grp1empty)) {//non-degenerate case
        Rcpp::NumericVector y1 = y_[grp1];
        Rcpp::NumericVector y2 = y_[!grp1];
        Rcpp::NumericVector x1 = beta_[grp1];
        Rcpp::NumericVector x2 = beta_[!grp1];
        b1 = sum(w1*y1)/sum(w1);
        b2 = sum(w2*(x2-y2))/sum(w2);
        xk = max(x1);
        xkp1 = min(x2);
        //compute unconstrained minimum for UB and LB
        a = b1 + b2;
        b = b1 - b2;
        //compute optimal UB and LB
        UB = std::max(std::min(xkp1, a), xk);
        LB = std::min(b, minval_);
        Rcpp::Rcout << "non-degenerate b1= " << b1 << " b2= " << b2 << " a= " << a << " b= " << b << " xk= " << xk << " xkp1= " << xkp1 << " UB= " << UB << " LB= " << LB << "\n";
    } else if (grp2empty) {//grp2 is empty, UB > xmax
        Rcpp::NumericVector w1 = weight_;
        Rcpp::NumericVector y1 = y_;
        b1 = sum(w1*y1)/sum(w1);
        b2 = 0;
        xk = xmax;
        xkp1 = std::numeric_limits<double>::infinity();
        //compute unconstrained minimum for UB and LB
        a = b1 + b2;
        b = b1 - b2;
        //compute optimal UB and LB
        if (b1 >= (minval_+xmax)/2.) {
            LB = minval_; //or any beta < minval_
            UB = 2*b1 - LB;
            Rcpp::Rcout << "case1";
        } else {
            UB = xmax; //or any beta > xmax
            LB = 2*b1 - UB;
            Rcpp::Rcout << "case2";
        }
        Rcpp::Rcout << "UB>xmax b1= " << b1 << " b2= " << b2 << " a= " << a << " b= " << b << " xk= " << xk << " xkp1= " << xkp1 << " UB= " << UB << " LB= " << LB << "\n";
    } else {//grp1 is empty, UB <= xmin
        if (!grp1empty) Rcpp::Rcout << "This should not have happened!\n";
        Rcpp::NumericVector w2 = weight_;
        Rcpp::NumericVector y2 = y_;
        Rcpp::NumericVector x2 = beta_;
        b1 = 0;
        b2 = sum(w2*(x2-y2))/sum(w2);
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
bounds_t BoundsOptimizer<ZeroOffset, Sign>::optimize_bounds(double val) const {
    //split data in two groups and determine constraint betas
    //Rcout << "  GET at " << val << " maxabsval_= " << maxabsval_ << std::endl;
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
    //Rcout << "  xk= " << xk << " xkp1= " << xkp1 << std::endl;
    //determine unconstrained UB
    double UB;
    //Rcout << "  grp1empty= " << grp1empty << " grp2empty= " << grp2empty << std::endl;
    if (grp1empty && grp2empty) {
        //UB larger than any beta
        UB=maxabsval_;
    } else if (grp1empty) {
        //UB larger than any positive beta
        Rcpp::NumericVector y2 = y_[grp2];
        Rcpp::NumericVector x2 = beta_[grp2];
        UB = -sum(w2*(x2-y2))/sum(w2);
    } else if (grp2empty) {
        Rcpp::NumericVector y1 = y_[grp1];
        Rcpp::NumericVector x1 = beta_[grp1];
        UB = sum(w1*(x1-y1))/sum(w1);
    } else {
        Rcpp::NumericVector y1 = y_[grp1];
        Rcpp::NumericVector x1 = beta_[grp1];
        double sw1 = sum(w1);
        double b1 = sum(w1*(x1-y1))/sw1;
        Rcpp::NumericVector y2 = y_[grp2];
        Rcpp::NumericVector x2 = beta_[grp2];
        double sw2 = sum(w2);
        double b2 = sum(w2*(x2-y2))/sw2;
        UB = (sw1*b1-sw2*b2)/(sw1+sw2);
    }
    //apply constraint
    UB = std::min(std::max(UB,xk),xkp1);
    return bounds_t{-UB,UB};
}

