#include <Rcpp.h>
#include <utility> //pair
#include <tuple> //tie

#include "Traits.hpp"
#include "ScoreAssembler.hpp"
#include "Likelihoods.hpp"
#include "DOFComputer.hpp"
#include "BinnedData.hpp"
#include "ScorePreparator.hpp"

template<typename Calculation, typename Score, typename GaussianEstimator>
ScoreComputer<Calculation,Score,GaussianEstimator>::ScoreComputer(double tol_val, const binned_t& data,
                                                GaussianEstimator& gauss, double lambda2)
 : ScorePreparator<Score,GaussianEstimator>(gauss, data, lambda2),
   Likelihood<Calculation>(data, ScorePreparator<Score,GaussianEstimator>::get_likelihood_var()),
   ScoreAssembler<Score>(ScorePreparator<Score,GaussianEstimator>::get_assembler_var()),
   DOFComputer(tol_val, data) {}


template<typename Calculation, typename Score, typename GaussianEstimator>
Rcpp::NumericVector ScoreComputer<Calculation,Score,GaussianEstimator>::evaluate(double LB, double UB) const {
    //compute dof and chi square
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    const Rcpp::NumericVector chisq = Likelihood<Calculation>::get_chi_square(LB, UB);
    const int dof = get_dof(LB, UB);
    
    //assemble it into a score
    double score, score_sd;
    std::tie(score, score_sd) = ScoreAssembler<Score>::assemble(chisq, dof);
    
    /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
     << " CV= " << score  << " dof= " << dof
     << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
    return Rcpp::NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                 _["BIC"]=score, _["BIC.sd"]=score_sd, _["dof"]=dof,
                                 _["UB"]=UB, _["LB"]=LB); // do not report score name for now
}

template<typename Calculation, typename Score, typename GaussianEstimator>
Rcpp::NumericVector ScoreComputer<Calculation,Score,GaussianEstimator>::invalidate(double LB, double UB) const {
    //compute dof and chi square
    const double lambda1 = (UB-LB)/2;
    const double eCprime = (UB+LB)/2.;
    const double score = std::numeric_limits<double>::max();
    const double score_sd = -1;
    /*Rcout << " OBJ " << msg << " ok lambda1= " << lambda1 << " eCprime= 0"
     << " CV= " << score  << " dof= NA"
     << " UB= " << lambda1 << " LB= " << -lambda1 << std::endl;*/
    return Rcpp::NumericVector::create(_["eCprime"]=eCprime, _["lambda1"]=lambda1,
                                       _["BIC"]=score, _["BIC.sd"]=score_sd, _["dof"]=NumericVector::get_na(),
                                       _["UB"]=UB, _["LB"]=LB); // do not report score name for now
}

