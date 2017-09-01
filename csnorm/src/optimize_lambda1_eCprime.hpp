#ifndef OPTIMIZE_LAMBDA1_ECPRIME_HPP
#define OPTIMIZE_LAMBDA1_ECPRIME_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <tuple> //tie

#include "Traits.hpp"
#include "ScoreComputer.hpp"
#include "BoundsComputer.hpp"
#include "Likelihoods.hpp"

class SignalData;

//objective functor to find lambda1 and eCprime assuming the signal is positive, using CV or BIC
template<typename Score>
class obj_lambda1_eCprime : private BoundsComputer<EstimatedOffset,PositiveSign>,
                             private ScoreComputer<Signal,Score> {
public:
    typedef typename ScoreComputer<Signal,Score>::likelihood_var_t likelihood_var_t;
    typedef typename ScoreComputer<Signal,Score>::assembler_var_t assembler_var_t;

    obj_lambda1_eCprime(double tol_val,
                        bool constrained, const SignalData& data, NumericVector forbidden_vals,
                        double lambda2, const likelihood_var_t& likelihood_var,
                        const assembler_var_t& assembler_var) :
       BoundsComputer<EstimatedOffset,PositiveSign>(data, min(data.get_value())),
       ScoreComputer<Signal,Score>(tol_val, data,  likelihood_var, assembler_var),
       tol_val_(tol_val), lambda2_(lambda2), constrained_(constrained), forbidden_vals_(forbidden_vals) {}
    
    NumericVector get(double val, std::string msg = "") const {
        //optimize bounds
        double UB, LB;
        std::tie(LB, UB) = BoundsComputer<EstimatedOffset,PositiveSign>::optimize_bounds(val);
        //check if forbidden. TODO: encapsulate
        double eCprime = (UB+LB)/2;
        double lambda1 = (UB-LB)/2;
        /*Rcout << "EVAL at val= " << val << " LB= " << LB << " UB= " << UB << " b1= " << b1 << " b2= " << b2
         << " xmin= " << xmin << " xk= " << xk << " xkp1= " << xkp1 << " a= " << a << " b= " << b << std::endl; */
        //check if solution is feasible
        if ( UB<LB || (constrained_ && is_true(any( (forbidden_vals_>UB+tol_val_/2) | (forbidden_vals_<LB-tol_val_/2) )) ) ) {
            if (!msg.empty()) Rcout << " OBJ " << msg << " forbidden lambda2= " << lambda2_ << " lambda1= " << lambda1
                << " eCprime= " << eCprime << " CV= Inf dof= NA"
                << " UB= " << UB  << " LB= " << LB << std::endl;
            return NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                         _["BIC"]=std::numeric_limits<double>::max(), _["BIC.sd"]=0,
                                         _["dof"]=NumericVector::get_na(), _["UB"]=UB, _["LB"]=LB);
        }
        //return score
        return ScoreComputer<Signal,Score>::evaluate(LB, UB);

    }

private:
    double tol_val_, lambda2_;
    bool constrained_;
    NumericVector forbidden_vals_;
};

NumericVector cpp_optimize_lambda1_eCprime(const DataFrame mat, int nbins,
        double tol_val, bool constrained,
        double lambda1_min, int refine_num, double lambda2);


#endif

