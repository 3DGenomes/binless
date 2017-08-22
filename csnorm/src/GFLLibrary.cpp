#include <Rcpp.h>
#include <vector>

#include "GFLLibrary.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm


void GFLLibrary::setUp(double alpha, double inflate, int ninner, double converge, double clamp) {
    inflate_ = inflate;
    ninner_ = ninner;
    converge_ = converge;
    clamp_ = clamp;
    
    alpha_ = alpha;
    counter_ = 0;
}

void GFLLibrary::optimize(const std::vector<double>& y, const std::vector<double>& beta_init,
                          const std::vector<double>& w, double lambda2) {
    //setup initial values
    beta_ = beta_init;
    std::vector<double> u(tsz_,0); //residuals set to zero
    std::vector<double> z;
    z.reserve(tsz_);
    for (int i=0; i<trails_.size(); ++i) {
        z.push_back(beta_[trails_[i]]);   //z set to beta values along trails
    }
    
    //perform optimization on the C side
    double* py = const_cast<double*>(&y[0]);
    double* pw = const_cast<double*>(&w[0]);
    counter_ += graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, &trails_[0], &breakpoints_[0],
                                               lambda2, &alpha_, inflate_, ninner_, converge_,
                                               &beta_[0], &z[0], &u[0]);
    
    //clamp values at +- clamp_ if needed
    for (std::vector<double>::iterator it = beta_.begin(); it != beta_.end(); ++it) {
        *it = std::min(clamp_, std::max(-clamp_, *it));
    }
    
}

std::vector<std::vector<int> > GFLLibrary::triangle_grid_chain(unsigned nrows) const {
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

