#ifndef OPTIMIZE_LAMBDA1_HPP
#define OPTIMIZE_LAMBDA1_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <tuple> // tie

#include "Traits.hpp"
#include "ScoreComputer.hpp"
#include "BoundsComputer.hpp"
#include "Traits.hpp"
#include "Likelihoods.hpp"

class SignalData;

//objective functor to find lambda1 assuming eCprime=0, using CV or BIC, signal case
template<typename Score>
class obj_lambda1 : private BoundsComputer<ZeroOffset,AnySign>,
                    private ScoreComputer<Signal,Score> {
public:
    obj_lambda1(double minUB, double tol_val,
                const SignalData& data, NumericVector forbidden_vals,
                const typename ScoreComputer<Signal,Score>::var_t& score_specific) :
      BoundsComputer<ZeroOffset,AnySign>(data, minUB),
      ScoreComputer<Signal,Score>(tol_val, data, score_specific),
      minUB_(minUB), tol_val_(tol_val), forbidden_vals_(forbidden_vals) {}
    
    NumericVector get(double val, std::string msg = "") const {
        //optimize bounds
        double UB,LB;
        std::tie(LB, UB) = BoundsComputer<ZeroOffset,AnySign>::optimize_bounds(val);
        //check if forbidden. TODO: encapsulate
        double lambda1=UB;
        if ( is_true(any(abs(forbidden_vals_)>lambda1+tol_val_/2)) || (UB < minUB_) ) {
            if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda1= " << lambda1
                << " eCprime= 0 CV= Inf dof= NA"
                << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;
            return NumericVector::create(_["eCprime"]=0, _["lambda1"]=lambda1,
                                         _["BIC"]=std::numeric_limits<double>::max(), _["BIC.sd"]=0,
                                         _["dof"]=NumericVector::get_na(), _["UB"]=lambda1, _["LB"]=-lambda1);
        }
        //return score
        return ScoreComputer<Signal,Score>::evaluate(LB, UB);
    }
    
private:
    double minUB_, tol_val_;
    NumericVector forbidden_vals_;
};


NumericVector cpp_optimize_lambda1(const DataFrame mat, int nbins,
                                   double tol_val, bool positive,
                                   double lambda1_min, int refine_num);


#endif

