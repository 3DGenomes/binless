#ifndef COMPUTE_CV_HPP
#define COMPUTE_CV_HPP


struct compute_CV_diff {
    compute_CV_diff(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
                    const NumericVector& weight_ref, const NumericVector& valuehat_ref,
                    const IntegerVector& patchno, const IntegerVector& cv_grp) :
    tol_val_(tol_val), value_(value), weight_(weight), valuehat_(valuehat), weight_ref_(weight_ref),
    valuehat_ref_(valuehat_ref), patchno_(patchno), cv_grp_(cv_grp) {}
    NumericVector evaluate(double LB, double UB) const;
private:
    double tol_val_;
    NumericVector value_, weight_, valuehat_, weight_ref_, valuehat_ref_;
    IntegerVector patchno_, cv_grp_;
};

struct compute_CV_signal {
    compute_CV_signal(double tol_val, const NumericVector& value, const NumericVector& weight, const NumericVector& valuehat,
                      const IntegerVector& patchno, const IntegerVector& cv_grp) :
    tol_val_(tol_val), value_(value), weight_(weight), valuehat_(valuehat), patchno_(patchno), cv_grp_(cv_grp) {}
    NumericVector evaluate(double LB, double UB) const;
private:
    double tol_val_;
    NumericVector value_, weight_, valuehat_;
    IntegerVector patchno_, cv_grp_;
};

#endif