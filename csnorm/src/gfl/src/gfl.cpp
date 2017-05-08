#include <Rcpp.h>
using namespace Rcpp;
#include "graph_fl.c"
#include "tf_dp.c"
#include "utils.c"

RcppExport SEXP weighted_graphfl(Rcpp::NumericVector y_i, Rcpp::NumericVector
        w_i, int ntrails, Rcpp::NumericVector trails_i, Rcpp::NumericVector breakpoints_i,
        double lam, double alpha, double inflate, int maxsteps, double converge)
{
  size_t N = y_i.size();
  std::vector<double> y_r = Rcpp::as<std::vector<double> >(y_i);
  std::vector<double> w_r = Rcpp::as<std::vector<double> >(w_i);
  std::vector<int> trails_r = Rcpp::as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = Rcpp::as<std::vector<int> >(breakpoints_i);
  std::vector<double> beta_r(N);
  std::vector<double> z_r(breakpoints_r[ntrails-1], 0);
  std::vector<double> u_r(breakpoints_r[ntrails-1], 0);
  
  int res;
  res = graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails, &trails_r[0], &breakpoints_r[0],
                                     lam, alpha, inflate, maxsteps, converge,
                                     &beta_r[0], &z_r[0], &u_r[0]);
  
  return Rcpp::wrap(beta_r);
}

RCPP_MODULE(gfl){
  using namespace Rcpp ;
 
  function("weighted_graphfl" , &weighted_graphfl  , "documentation for weighted_graphfl ");
} 

