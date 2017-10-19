#ifndef LIKELIHOODS_HPP
#define LIKELIHOODS_HPP

#include <Rcpp.h>

#include "BinnedData.hpp"
#include "Traits.hpp"
#include "util.hpp"
#include "Settings.hpp"

//Likelihood takes the data and unthresholded beta estimate, and returns
//the squared deviations on each matrix bin, computed on the
//templated data type (Signal or Difference)
//(commented out because already declared in Settings)
//template<typename Calculation> class Likelihood {};

//Signal
template<> class Likelihood<Signal> {
public:
    typedef Rcpp::NumericVector var_t;
    
    Likelihood(const BinnedData<Signal>& data, const var_t& beta_phi)
      : binned_(data), beta_phi_(beta_phi) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const {
        //soft-threshold to get phi
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> beta_phi_r = Rcpp::as<std::vector<double> >(beta_phi_);
        Rcpp::NumericVector phi = Rcpp::wrap(soft_threshold(beta_phi_r, eCprime, lambda1));
        
        //cap the weights
        Rcpp::NumericVector wt = binned_.get_weight();
        const double wmax = get_maximum_admissible_weight(wt, Settings<Likelihood<Signal> >::get_weights_max_percentile());
        Rcpp::Rcout << "maximum admissible weight: " << wmax << "\n";
        wt = Rcpp::pmin(wt, wmax);
        
        //compute the score and return
        const Rcpp::NumericVector chisq = wt * SQUARE(binned_.get_phihat() - (phi + eCprime));
        
        /*Rcpp::Rcout << "UB= " << UB << " LB= " << LB << " eCprime= " << eCprime << " lambda1= " << lambda1
            << " max(beta)= " << (double)(Rcpp::max(beta_phi_))
            << " min(beta)= " << (double)(Rcpp::min(beta_phi_))
            << " max(phi)= " << (double)(Rcpp::max(phi))
            << " min(phi)= " << (double)(Rcpp::min(phi))
            << " sum(chisq)= " << (double)(Rcpp::sum(chisq)) << "\n";
            Rcpp::Rcout << "beta_phi:\n";
            for (auto i : Rcpp::as<std::vector<double> >(beta_phi_)) Rcpp::Rcout << i << "\n";
            Rcpp::Rcout << "done beta_phi\n";*/
        return chisq;
    }
    
private:
    const BinnedData<Signal>& binned_;
    const Rcpp::NumericVector beta_phi_;
};

//Difference
template<> class Likelihood<Difference> {
public:
    typedef Rcpp::NumericVector var_t;
    
    Likelihood(const BinnedData<Difference>& data, const var_t& beta_delta)
      : binned_(data), beta_delta_(beta_delta) {}
    
    Rcpp::NumericVector get_chi_square(double LB, double UB) const {
        //soft-threshold to get delta and phi_ref
        const double lambda1 = (UB-LB)/2;
        const double eCprime = (UB+LB)/2.;
        std::vector<double> beta_delta_r = Rcpp::as<std::vector<double> >(beta_delta_);
        std::vector<double> delta_r = soft_threshold(beta_delta_r, eCprime, lambda1);
        Rcpp::NumericVector delta = Rcpp::wrap(delta_r);
        Rcpp::NumericVector phi_ref = compute_phi_ref(binned_, delta);
        
        //cap the weights
        Rcpp::NumericVector wt = binned_.get_weight();
        const double wmax = get_maximum_admissible_weight(wt, Settings<Likelihood<Signal> >::get_weights_max_percentile());
        wt = Rcpp::pmin(wt, wmax);
        Rcpp::NumericVector wt_ref = binned_.get_weight_ref();
        const double wmax_ref = get_maximum_admissible_weight(wt_ref, Settings<Likelihood<Signal> >::get_weights_max_percentile());
        wt_ref = Rcpp::pmin(wt_ref, wmax_ref);
        Rcpp::Rcout << "maximum admissible weights: " << wmax << " " << wmax_ref << "\n";
        
        const Rcpp::NumericVector chisq = wt * SQUARE(binned_.get_phihat() - (delta + eCprime + phi_ref))
                                        + wt_ref * SQUARE(binned_.get_phihat_ref() - phi_ref);
        return chisq;
    }
    
private:
    const BinnedData<Difference>& binned_;
    const Rcpp::NumericVector beta_delta_;
};


#endif
