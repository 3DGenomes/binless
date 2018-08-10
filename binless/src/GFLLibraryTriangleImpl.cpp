#include <Rcpp.h>
#include <vector>

#include "GFLLibraryTriangleImpl.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm

std::vector<std::vector<int> > triangle_grid_chain(int nrows) {
    unsigned ntotal = nrows;
    ntotal = ntotal*(ntotal+1)/2-1;
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

void GFLLibraryTriangleImpl::store_trails(int nrows) {
    const std::vector<std::vector<int> > chains = triangle_grid_chain(nrows);
    //
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

void GFLLibraryTriangleImpl::reset() {
    //setup initial values for a cold start
    set_ninner(0);
    beta_impl_ = std::vector<double>(N_,0);
    z_ = std::vector<double>(tsz_,0);
    u_ = std::vector<double>(tsz_,0);
}

void GFLLibraryTriangleImpl::optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
    //perform optimization on the C side
    double* py = const_cast<double*>(&y[0]);
    double* pw = const_cast<double*>(&w[0]);
    double alpha = get_alpha();
    int counter = graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, &trails_[0], &breakpoints_[0],
                                               lambda2, &alpha, get_inflate(), get_ninner_max(), converge,
                                               &beta_impl_[0], &z_[0], &u_[0]);
    set_alpha(alpha);
    set_beta(beta_impl_);
    set_ninner(get_ninner() + counter);
}

GFLLibraryTriangleImpl::GFLState_t GFLLibraryTriangleImpl::get_state() const {
    return Rcpp::List::create(_["z"]=z_, _["u"]=u_,  _["alpha"]=get_alpha(),  _["beta"]=beta_impl_,
                              _["counter"]=get_ninner(), _["maxdiag"]=0);
}

void GFLLibraryTriangleImpl::set_state(const GFLState_t& state) {
    if (state.containsElementNamed("u") && Rcpp::as<std::vector<double> >(state["u"]).size() == tsz_
          && state.containsElementNamed("maxdiag") && Rcpp::as<int>(state["maxdiag"])==0 ) {
        z_ = Rcpp::as<std::vector<double> >(state["z"]);
        u_ = Rcpp::as<std::vector<double> >(state["u"]);
        beta_impl_ = Rcpp::as<std::vector<double> >(state["beta"]);
        set_alpha(Rcpp::as<double>(state["alpha"]));
        set_ninner(Rcpp::as<unsigned>(state["counter"]));
    }
}


