#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <algorithm>

#include "fast_binless.hpp"
#include "GFLLibrary.hpp"
#include "FusedLassoGaussianEstimator.hpp"
 
std::vector<double> fast_compute_exposures(const FastSignalData& data) {
    //get residuals
    ResidualsPair z = get_poisson_residuals(data);
    //average by name
    std::vector<double> exposures(data.get_ndatasets(),0);
    std::vector<double> weightsums(data.get_ndatasets(),0);
    auto names = data.get_name();
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned name = names[i]-1; //offset by 1 for vector indexing
        exposures[name] += z.residuals[i]*z.weights[i];
        weightsums[name] += z.weights[i];
    }
    //add current estimates
    auto expo_ori = data.get_exposures();
    for (unsigned i=0; i<data.get_ndatasets(); ++i) {
        exposures[i] = exposures[i]/weightsums[i] + expo_ori[i];
    }
    return exposures;
}

std::vector<double> fast_compute_log_biases(const FastSignalData& data) {
    //get residuals
    ResidualsPair z = get_poisson_residuals(data);
    //sum them along the rows/columns
    auto dbin1 = data.get_bin1();
    auto dbin2 = data.get_bin2();
    std::vector<double> log_biases(data.get_nbins(),0);
    std::vector<double> weightsums(data.get_nbins(),0);
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = dbin2[i]-1;
        double w = z.weights[i];
        double b = z.residuals[i]*w;
        log_biases[bin1] += b;
        weightsums[bin1] += w;
        if (bin1!=bin2) {
            log_biases[bin2] += b;
            weightsums[bin2] += w;
        }
    }
    //add current bias and compute weighted average
    double avg=0;
    double wsum=0;
    auto ori_log_biases = data.get_log_biases();
    for (unsigned i=0; i<data.get_nbins(); ++i) {
        log_biases[i] = (weightsums[i]>0) ? log_biases[i]/weightsums[i] : 0;
        log_biases[i] += ori_log_biases[i];
        avg += log_biases[i]*weightsums[i];
        wsum += weightsums[i];
    }
    avg = avg / wsum;
    //no smoothing, subtract average and return
    for (unsigned i=0; i<data.get_nbins(); ++i) {
        log_biases[i] -= avg;
    }
    return log_biases;
}


std::vector<double> fast_compute_log_decay(const FastSignalData& data) {
    //get residuals
    ResidualsPair z = get_poisson_residuals(data);
    //sum them along the diagonals
    auto dbin1 = data.get_bin1();
    auto dbin2 = data.get_bin2();
    std::vector<double> log_decay(data.get_nbins(),0);
    std::vector<double> weightsums(data.get_nbins(),0);
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = dbin2[i]-1;
        log_decay[bin2-bin1] += z.residuals[i]*z.weights[i];
        weightsums[bin2-bin1] += z.weights[i];
    }
    //add current bias and compute weighted average
    double avg=0;
    double wsum=0;
    auto dlog_decay = data.get_log_decay();
    for (unsigned i=0; i<data.get_nbins(); ++i) {
        log_decay[i] = (weightsums[i]>0) ? log_decay[i]/weightsums[i] : 0;
        log_decay[i] += dlog_decay[i];
        avg += log_decay[i]*weightsums[i];
        wsum += weightsums[i];
    }
    avg = avg / wsum;
    //no smoothing, subtract average and return
    for (unsigned i=0; i<data.get_nbins(); ++i) {
        log_decay[i] -= avg;
    }
    return log_decay;
}

double fast_precision(const std::vector<double>& weights, const std::vector<double>& weights_old) {
    double delta = std::abs(weights[0]-weights_old[0]);
    double maxval = weights[0];
    double minval = weights[0];
    const unsigned N = weights.size();
    for (unsigned i=1; i < N; ++i) {
        delta = std::max(std::abs(weights[i]-weights_old[i]), delta);
        minval = std::min(weights[i],minval);
        maxval = std::max(weights[i],maxval);
    }
    return delta/(maxval-minval);
}

std::vector<double> fast_remove_signal_degeneracy(const FastSignalData& data) {
    auto dbin1 = data.get_bin1();
    auto dbin2 = data.get_bin2();
    auto dname = data.get_name();
    std::vector<double> log_signal = data.get_log_signal();
    double max_signal = *std::max_element(log_signal.begin(),log_signal.end());
    //get minimum signal per row and per counter diagonal
    std::vector<std::vector<double> > min_per_row(data.get_nbins(),std::vector<double>(data.get_ndatasets(),max_signal));
    std::vector<std::vector<double> > min_per_diag(data.get_nbins(),std::vector<double>(data.get_ndatasets(),max_signal));
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned name = dname[i]-1;
        unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = dbin2[i]-1;
        min_per_diag[bin2-bin1][name] = std::min(min_per_diag[bin2-bin1][name], log_signal[i]);
        min_per_row[bin1][name] = std::min(min_per_row[bin1][name], log_signal[i]);
        if (bin1!=bin2) {
            min_per_row[bin2][name] = std::min(min_per_row[bin2][name], log_signal[i]);
        }
    }
    //remove from each signal row and counter diagonal
    double max_adjust;
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned name = dname[i]-1;
        unsigned bin1 = dbin1[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = dbin2[i]-1;
        double adjust = std::max(std::max(min_per_row[bin1][name],min_per_row[bin2][name]),min_per_diag[bin2-bin1][name]);
        log_signal[i] = std::max(log_signal[i]-adjust,0.);
        max_adjust = (i==0) ? adjust : std::max(adjust,max_adjust);
    }
    Rcpp::Rcout << " max adjustment for constraint: " << max_adjust << "\n";
    return log_signal;
}

std::vector<double> fast_shift_signal(const FastSignalData& data) {
    auto dname = data.get_name();
    std::vector<double> log_signal = data.get_log_signal();
    double max_signal = *std::max_element(log_signal.begin(),log_signal.end());
    //get minimum signal per dataset
    std::vector<double> min_per_dset(data.get_ndatasets(),max_signal);
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned name = dname[i]-1;
        min_per_dset[name] = std::min(min_per_dset[name], log_signal[i]);
    }
    //shift signals
    for (unsigned i=0; i<data.get_N(); ++i) {
        unsigned name = dname[i]-1;
        log_signal[i] = std::max(log_signal[i]-min_per_dset[name],0.);
    }
    return log_signal;
}

List fast_binless(const DataFrame obs, unsigned nbins, unsigned ngibbs, double lam2, double tol_val) {
    //initialize return values, exposures and fused lasso optimizer
    Rcpp::Rcout << "init\n";
    FastSignalData out(obs, nbins);
    out.set_exposures(fast_compute_exposures(out));
    const double converge = tol_val/20.;
    std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(out.get_ndatasets(),
                                                               FusedLassoGaussianEstimator<GFLLibrary>(nbins, converge));
    std::vector<DataFrame> diagnostics;
    for (unsigned step=1; step <= ngibbs; ++step) {
        Rcpp::Rcout << "step " << step << "\n";
        //compute biases
        Rcpp::Rcout << " biases\n";
        auto biases = fast_compute_log_biases(out);
        out.set_log_biases(biases);
        //compute decay
        Rcpp::Rcout << " decay\n";
        auto decay = fast_compute_log_decay(out);
        out.set_log_decay(decay);
        //compute signal
        Rcpp::Rcout << " signal\n";
        auto old_weights = out.get_signal_weights();
        auto signal = fast_compute_signal(out, flos, lam2);
        out.set_log_signal(signal.beta);
        out.set_signal_phihat(signal.phihat);
        out.set_signal_weights(signal.weights);
        //compute precision and convergence
        double precision = fast_precision(out.get_signal_weights(),old_weights);
        Rcpp::Rcout << " reached relative precision " << precision << "\n";
        bool converged = precision < tol_val;
        if (converged || step == ngibbs) {
            auto adjust = fast_shift_signal(out);
            out.set_log_signal(adjust);
        } else {
            auto adjust = fast_remove_signal_degeneracy(out);
            out.set_log_signal(adjust);
        }
        //compute exposures
        Rcpp::Rcout << " exposures\n";
        auto exposures = fast_compute_exposures(out);
        out.set_exposures(exposures);
        diagnostics.push_back(out.get_as_dataframe());
        if (converged) {
            Rcpp::Rcout << "converged or reached last step\n";
            break;
        }
    }
    Rcpp::Rcout << "done\n";
    //finalize and return
    return Rcpp::List::create(_["mat"]=out.get_as_dataframe(), _["log_biases"]=out.get_log_biases(),
                              _["log_decay"]=out.get_log_decay(), _["exposures"]=out.get_exposures(),
                              _["diagnostics"]=diagnostics, _["nbins"]=nbins);
}
 
Rcpp::DataFrame fast_binless_difference(const List obs, double lam2, double tol_val, unsigned ref) {
    //read normalization data
    const unsigned nbins = obs["nbins"];
    const Rcpp::DataFrame mat = Rcpp::as<Rcpp::DataFrame>(obs["mat"]);
    FastDifferenceData out(mat, nbins, ref);
    auto signal = Rcpp::as<std::vector<double> >(mat["signal"]);
    out.set_log_signal(signal); //fills-in phi_ref and delta
    auto log_decay = Rcpp::as<std::vector<double> >(obs["log_decay"]);
    out.set_log_decay(log_decay);
    auto log_biases = Rcpp::as<std::vector<double> >(obs["log_biases"]);
    out.set_log_biases(log_biases);
    auto exposures = Rcpp::as<std::vector<double> >(obs["exposures"]);
    out.set_exposures(exposures);
    const double converge = tol_val/20.;
    std::vector<FusedLassoGaussianEstimator<GFLLibrary> > flos(out.get_ndatasets(),
                                                               FusedLassoGaussianEstimator<GFLLibrary>(nbins, converge));
    //compute differences
    auto old_weights = out.get_difference_weights();
    auto diff = fast_compute_difference(out, flos, lam2, ref);
    out.set_log_difference(diff.delta);
    out.set_deltahat(diff.deltahat);
    out.set_difference_weights(diff.weights);
    out.set_phi_ref(diff.phi_ref);
    //finalize and return
    return out.get_as_dataframe();
}

