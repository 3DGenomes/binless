#include <assert.h>
#include <Rcpp.h>
using namespace Rcpp;
#include "graph_fl.c"
#include "tf_dp.c"
#include "utils.c"
#include "gfl_c.c"

RcppExport SEXP weighted_graphfl(Rcpp::NumericVector y_i, Rcpp::NumericVector
        w_i, int ntrails, Rcpp::NumericVector trails_i, Rcpp::NumericVector breakpoints_i,
        double lam, double alpha, double inflate, int ninner, double converge)
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
                                       lam, &alpha, inflate, ninner, converge,
                                       &beta_r[0], &z_r[0], &u_r[0]);
  
  return Rcpp::wrap(beta_r);
}


// [[Rcpp::export]]
DataFrame cts_to_signal_mat(const DataFrame cts, int nbins, double dispersion, std::vector<double>& phi, int diag_rm)
{
  //inputs
  int N = cts.nrows();
  std::vector<int> cts_bin1 = as<std::vector<int> >(cts["bin1"]);
  std::vector<int> cts_bin2 = as<std::vector<int> >(cts["bin2"]);
  std::vector<double> count = as<std::vector<double> >(cts["count"]);
  std::vector<double> lmu_nosig = as<std::vector<double> >(cts["lmu.nosig"]);
  std::vector<double> weight = as<std::vector<double> >(cts["weight"]);
  
  //outputs
  int nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<double> phihat(nbetas, 0); //vectorized form
  std::vector<double> phihat_var(nbetas, 0);
  std::vector<double> ncounts(nbetas, 0);
  std::vector<int> bin1(nbetas, 0);
  std::vector<int> bin2(nbetas, 0);
  
  cts_to_signal_mat_core(N, &cts_bin1[0], &cts_bin2[0], &count[0], &lmu_nosig[0], &weight[0], nbins, dispersion, &phi[0],
                  &phihat[0], &phihat_var[0], &ncounts[0], &bin1[0], &bin2[0], diag_rm);
  
  IntegerVector bin1_i, bin2_i;
  NumericVector phihat_i, phihat_var_i, ncounts_i, weight_i, didx_i;
  bin1_i = wrap(bin1);
  bin2_i = wrap(bin2);
  bin1_i.attr("levels") = as<IntegerVector>(cts["bin1"]).attr("levels");
  bin2_i.attr("levels") = as<IntegerVector>(cts["bin2"]).attr("levels");
  bin1_i.attr("class") = CharacterVector::create("ordered", "factor");
  bin2_i.attr("class") = CharacterVector::create("ordered", "factor");
  phihat_i = wrap(phihat);
  phihat_var_i = wrap(phihat_var);
  ncounts_i = wrap(ncounts);
  weight_i = 1/phihat_var_i;
  didx_i = bin2_i-bin1_i;
  
  
  return DataFrame::create(_["bin1"]=bin1_i, _["bin2"]=bin2_i, _["phihat"]=phihat_i,
                           _["phihat.var"]=phihat_var_i, _["ncounts"]=ncounts_i, _["weight"]=weight_i,
                           _["diag.idx"]=didx_i);
}

// [[Rcpp::export]]
List wgfl_signal_perf_warm(const DataFrame cts, double dispersion, int nouter, int nbins,
                    int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                    double lam,  double alpha, double inflate, int ninner, double converge,
                    int diag_rm, NumericVector z_i, NumericVector u_i, NumericVector phi_i)
{
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<int> trails_r = as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
  std::vector<double> z_r = as<std::vector<double> >(z_i);
  std::vector<double> u_r = as<std::vector<double> >(u_i);
  std::vector<double> phi_r = as<std::vector<double> >(phi_i);
  std::vector<double> phi_old = phi_r;
  
  int step;
  int res=0;
  //printf(" Perf iteration: start with lam=%f alpha=%f phi[0]=%f z[0]=%f u[0]=%f\n",
  //       lam, alpha, phi_r[0], z_r[0], u_r[0]);
  for (step=1; step<=nouter; ++step) {
    const DataFrame mat = cts_to_signal_mat(cts, nbins, dispersion, phi_r, diag_rm);
    std::vector<double> y_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
    std::vector<double> w_r = Rcpp::as<std::vector<double> >(mat["weight"]);
    
    res += graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails, &trails_r[0], &breakpoints_r[0],
                                          lam, &alpha, inflate, ninner, converge,
                                          &phi_r[0], &z_r[0], &u_r[0]);
    double maxval = std::abs(phi_r[0]-phi_old[0]);
    for (int i=1; i<N; ++i) maxval = std::max(std::abs(phi_r[i]-phi_old[i]), maxval);
    //printf(" Iteration %d with alpha=%f reached maxval=%.5e after %d steps (phi[0]=%f z[0]=%f u[0]=%f)\n",
    //       step,alpha,maxval,res,phi_r[0],z_r[0], u_r[0]);
    if (maxval<converge) break;
    phi_old = phi_r;
  }
  //printf(" Perf iteration: end   with lam=%f alpha=%f phi[0]=%f z[0]=%f u[0]=%f nouter=%d ninner=%d\n",
  //       lam, alpha, phi_r[0], z_r[0], u_r[0], step, res);
  return List::create(_["phi"]=wrap(phi_r), _["alpha"]=wrap(alpha),
                      _["mat"]=cts_to_signal_mat(cts, nbins, dispersion, phi_r, diag_rm),
                      _["z"]=wrap(z_r), _["u"]=wrap(u_r), _["nouter"]=step, _["ninner"]=res);
}

// [[Rcpp::export]]
List wgfl_signal_perf(const DataFrame cts, double dispersion, int nouter, int nbins,
        int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
        double lam,  double alpha, double inflate, int ninner, double converge, int diag_rm)
{
    NumericVector z_i(breakpoints_i(ntrails-1));
    NumericVector u_i(breakpoints_i(ntrails-1));
    const int N = nbins*(nbins+1)/2; //size of fused lasso problem
    NumericVector phi_i(N);
    //printf("Fused lasso cold perf iteration with %d coefficients\n",phi.size());
    return wgfl_signal_perf_warm(cts, dispersion, nouter, nbins, ntrails, trails_i, breakpoints_i,
                                 lam, alpha, inflate, ninner, converge, diag_rm, z_i, u_i, phi_i);
}


RCPP_MODULE(gfl){
  using namespace Rcpp ;
 
  function("weighted_graphfl" , &weighted_graphfl  , "documentation for weighted_graphfl ");
  function("cts_to_signal_mat" , &cts_to_signal_mat  , "documentation for cts_to_signal_mat ");
  function("wgfl_signal_perf" , &wgfl_signal_perf  , "documentation for wgfl_signal_perf ");
  function("wgfl_signal_perf_warm" , &wgfl_signal_perf_warm  , "documentation for wgfl_signal_perf_warm ");
} 

