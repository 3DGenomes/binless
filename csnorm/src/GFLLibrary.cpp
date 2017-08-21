#include <Rcpp.h>
#include <vector>

#include "GFLLibrary.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm


void GFLLibrary::setUp(int ntrails, const NumericVector& trails, const NumericVector& breakpoints,
                  double alpha, double inflate, int ninner, double converge, double clamp) {
    ntrails_ = ntrails;
    trails_ = Rcpp::as<std::vector<int> >(trails);
    breakpoints_ = Rcpp::as<std::vector<int> >(breakpoints);
    tsz_ = trails_.size();
    inflate_ = inflate;
    ninner_ = ninner;
    converge_ = converge;
    clamp_ = clamp;
    
    alpha_ = alpha;
    counter_ = 0;
}

void GFLLibrary::optimize(const std::vector<double>& y, const std::vector<double>& beta_init,
                          const std::vector<double>& w, double lambda2) const {
    //setup initial values
    beta_ = beta_init;
    const unsigned tsz=trails_.size(); //beware that ntrails_ is not tsz
    std::vector<double> u(tsz,0); //residuals set to zero
    std::vector<double> z;
    z.reserve(tsz);
    for (int i=0; i<trails_.size(); ++i) {
        z.push_back(beta_[trails_[i]]);   //z set to beta values along trails
    }
    
    //perform optimization on the C side
    double* py = const_cast<double*>(&y[0]);
    double* pw = const_cast<double*>(&w[0]);
    int* ptr = const_cast<int*>(&trails_[0]);
    int* pbr = const_cast<int*>(&breakpoints_[0]);
    counter_ += graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, ptr, pbr,
                                               lambda2, &alpha_, inflate_, ninner_, converge_,
                                               &beta_[0], &z[0], &u[0]);
    
    //clamp values at +- clamp_ if needed
    for (std::vector<double>::iterator it = beta_.begin(); it != beta_.end(); ++it) {
        *it = std::min(clamp_, std::max(-clamp_, *it));
    }
    
}
