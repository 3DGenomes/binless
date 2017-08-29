#include <Rcpp.h>
using namespace Rcpp;

#include "base_objectives.hpp"
#include "util.hpp"

    
obj_lambda1_eCprime_base::bounds_t obj_lambda1_eCprime_base::optimize_bounds(double val) const {
    //split data in two groups
    LogicalVector grp1 = value_ <= val;
    NumericVector w1 = weight_[grp1];
    NumericVector w2 = weight_[!grp1];
    bool grp1only = w2.size() == 0 || sum(w2) == 0;
    bool grp2only = w1.size() == 0 || sum(w1) == 0;
    double b1,b2, xk, xkp1;
    double a,b,UB,LB;
    double xmin = min(value_), xmax = max(value_); //here, we assume xmin >= minval_
  
  //compute UB and LB
  if ((!grp1only) && (!grp2only)) {//non-degenerate case
      NumericVector betahat1 = valuehat_[grp1];
      NumericVector betahat2 = valuehat_[!grp1];
      NumericVector x1 = value_[grp1];
      NumericVector x2 = value_[!grp1];
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
  } else if (grp1only) {//grp2 is empty, UB > xmax
      NumericVector w1 = weight_;
      NumericVector betahat1 = valuehat_;
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
      if (!grp2only) Rcout << "This should not have happened!\n";
      NumericVector w2 = weight_;
      NumericVector betahat2 = valuehat_;
      NumericVector x2 = value_;
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
  return bounds_t{LB, UB};
}


obj_lambda1_base::bounds_t obj_lambda1_base::optimize_bounds(double val) const {
    //split data in two groups and determine constraint values
    //Rcout << "  GET at " << val << " maxabsval_= " << maxabsval_ << std::endl;
    LogicalVector grp1 = value_ > val, grp2 = value_ < -val;
    double xk,xkp1;
    if (val>maxabsval_) {
        xk=maxabsval_;
        xkp1=std::numeric_limits<double>::infinity();
    } else if (val<minabsval_) {
        xk=0;
        xkp1=minabsval_;
    } else {
        xk=max(as<NumericVector>(absval_[absval_<=val]));
        xkp1=min(as<NumericVector>(absval_[absval_>val]));
    }
    //Rcout << "  xk= " << xk << " xkp1= " << xkp1 << std::endl;
    //determine unconstrained UB
    double UB;
    bool grp1empty = is_true(all(!grp1)), grp2empty = is_true(all(!grp2));
    //Rcout << "  grp1empty= " << grp1empty << " grp2empty= " << grp2empty << std::endl;
    if (grp1empty && grp2empty) {
        //UB larger than any value
        UB=maxabsval_;
    } else if (grp1empty) {
        //UB larger than any positive value
        NumericVector w2 = weight_[grp2];
        NumericVector betahat2 = valuehat_[grp2];
        NumericVector x2 = value_[grp2];
        UB = -sum(w2*(x2-betahat2))/sum(w2);
    } else if (grp2empty) {
        NumericVector w1 = weight_[grp1];
        NumericVector betahat1 = valuehat_[grp1];
        NumericVector x1 = value_[grp1];
        UB = sum(w1*(x1-betahat1))/sum(w1);
    } else {
        NumericVector w1 = weight_[grp1];
        NumericVector betahat1 = valuehat_[grp1];
        NumericVector x1 = value_[grp1];
        double sw1 = sum(w1);
        double b1 = sum(w1*(x1-betahat1))/sw1;
        NumericVector w2 = weight_[grp2];
        NumericVector betahat2 = valuehat_[grp2];
        NumericVector x2 = value_[grp2];
        double sw2 = sum(w2);
        double b2 = sum(w2*(x2-betahat2))/sw2;
        UB = (sw1*b1-sw2*b2)/(sw1+sw2);
    }
    //apply constraint
    UB = std::min(std::max(UB,xk),xkp1);
    if (minUB_ <= xkp1 && minUB_ > xk) UB=std::max(minUB_,UB);
    //Rcout << "  UB= " << UB << " minUB= " << minUB_ << std::endl;
    return UB;
}

