#ifndef OPTIMIZE_LAMBDA1_DIFF_HPP
#define OPTIMIZE_LAMBDA1_DIFF_HPP

#include <Rcpp.h>
using namespace Rcpp;

#include "CalculationTraits.hpp"
#include "ScoreComputer.hpp"
#include "Bounds.hpp"

class DifferenceData;

//objective functor to find lambda1 assuming eCprime=0, using CV or BIC, difference case
template<typename Score>
class obj_lambda1_diff : private BoundsComputer<ZeroOffset,AnySign>,
                         private ScoreComputer<Difference,Score> {
public:
    obj_lambda1_diff(double minUB, double tol_val,
                        const DifferenceData& data, NumericVector forbidden_vals,
                        const typename Score::var_t& score_specific) :
    BoundsComputer<ZeroOffset,AnySign>(data, minUB),
    ScoreComputer<Difference,Score>(tol_val, data, score_specific),
    minUB_(minUB), tol_val_(tol_val), forbidden_vals_(forbidden_vals) {}
    
    
    NumericVector get(double val, std::string msg = "") const {
        //optimize bounds
        double UB = BoundsComputer<ZeroOffset,AnySign>::optimize_bounds(val);
        //check if forbidden. TODO: encapsulate
        double lambda1=UB;
        if ( is_true(any(abs(forbidden_vals_)>lambda1+tol_val_/2)) || (UB < minUB_) ) {
            if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1
                << " eCprime= 0 CV= Inf dof= NA"
                << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
            return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                         _["BIC"]=std::numeric_limits<double>::max(), _["BIC.sd"]=0,
                                         _["dof"]=NumericVector::get_na(),
                                         _["UB"]=lambda1, _["LB"]=-lambda1);
        }
        //return score
        return ScoreComputer<Difference,Score>::evaluate(-UB, UB);
    }
    
    double minUB_, tol_val_;
    NumericVector forbidden_vals_;
};




NumericVector cpp_optimize_lambda1_diff(const DataFrame mat, int nbins,
                                   double tol_val, double lambda1_min, int refine_num);


#endif

