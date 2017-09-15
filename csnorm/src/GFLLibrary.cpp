#include <Rcpp.h>
#include <vector>

#include "GFLLibrary.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm

std::vector<std::vector<int> > triangle_grid_chain(unsigned nrows) {
    int ntotal = nrows*(nrows+1)/2-1;
    std::vector<std::vector<int> > chains;
    int l = nrows;
    std::vector<int> current(1,0);
    //rows of consecutive numbers
    for (int i=1; i<=ntotal; ++i) {
        if (current.size()==l) {
            chains.push_back(current);
            current = std::vector<int>(1,i);
            l--;
        } else {
            current.push_back(i);
        }
    }
    //columns with Ui+1 = Ui + (N-i) with U1 from 2 to nrow
    for (int U1=2; U1<=nrows; ++U1) {
        int Ui=U1;
        current = std::vector<int>(1,Ui-1);
        for (int i=1; i<U1; ++i) {
            int Uip1 = Ui + nrows - i;
            current.push_back(Uip1-1);
            Ui=Uip1;
        }
        chains.push_back(current);
    }
    return(chains);
}

void GFLLibrary::store_trails(unsigned nrows) {
    const std::vector<std::vector<int> > chains = triangle_grid_chain(nrows);
    trails_.clear();
    breakpoints_.clear();
    for (std::vector<std::vector<int> >::const_iterator it = chains.begin() ;
         it != chains.end(); ++it) {
        if (trails_.size()>0) breakpoints_.push_back(trails_.size());
        trails_.insert(trails_.end(), it->begin(), it->end());
    }
    if (trails_.size()>0) breakpoints_.push_back(trails_.size());
    ntrails_ = breakpoints_.size();
    tsz_ = trails_.size();
}

void GFLLibrary::reset() {
    //setup initial values for a cold start
    counter_ = 0;
    beta_ = std::vector<double>(N_,0);
    z_ = std::vector<double>(tsz_,0);
    u_ = std::vector<double>(tsz_,0);
}

void GFLLibrary::optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2) {
    //perform optimization on the C side
    double* py = const_cast<double*>(&y[0]);
    double* pw = const_cast<double*>(&w[0]);
    int counter = graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, &trails_[0], &breakpoints_[0],
                                               lambda2, &alpha_, inflate_, ninner_, converge_,
                                               &beta_[0], &z_[0], &u_[0]);
    //Rcpp::Rcout << "GFLLibrary: " << counter << " steps\n";
    counter_ += counter;
}

GFLLibrary::GFLState_t GFLLibrary::get_state() const {
    return Rcpp::List::create(_["z"]=z_, _["u"]=u_,  _["alpha"]=alpha_,  _["beta"]=beta_,  _["counter"]=counter_);
}

void GFLLibrary::set_state(const GFLState_t& state) {
    if (state.containsElementNamed("u") && Rcpp::as<std::vector<double> >(state["u"]).size() == tsz_) {
        z_ = Rcpp::as<std::vector<double> >(state["z"]);
        u_ = Rcpp::as<std::vector<double> >(state["u"]);
        beta_ = Rcpp::as<std::vector<double> >(state["beta"]);
        alpha_ = Rcpp::as<double>(state["alpha"]);
        counter_ = Rcpp::as<unsigned>(state["counter"]);
    }
}


