#include <assert.h>
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

void cts_to_mat_core(int N, int* cts_bin1, int* cts_bin2, double* count, double* lmu_nosig, double* weight,
                          int nbins, double dispersion, double* phi, double* phihat, double* phihat_var,
                          double* ncounts, int* bin1, int* bin2)
{
  //walk through cts
  int nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
  int i;
  for (i = 0; i < N; ++i) {
    //pos(b1,b2) starts at 0 and is the index in the 2D triangle grid
    //of the cell at coordinates (b1,b2), where b1 and b2 start at 1 and b2>=b1
    int b1 = cts_bin1[i];
    int b2 = cts_bin2[i];
    int pos = (b1-1)*(nbins+1) - (b1-1)*b1/2 + b2 - b1;
    double mu = std::exp( lmu_nosig[i] + phi[pos] );
    double z = count[i]/mu-1;
    double var = 1/mu + 1/dispersion;
    double w2v = weight[i]/(2*var);
    ncounts[pos] += weight[i];
    phihat_var[pos] += w2v;
    phihat[pos] += (z+phi[pos])*w2v;
  }
  
  //finish mat
  int b1 = 1;
  int b2 = 1;
  for (i = 0; i < nbetas; ++i) {
    phihat_var[i] = 1/phihat_var[i];
    phihat[i] = phihat[i]*phihat_var[i];
    bin1[i] = b1;
    bin2[i] = b2++;
    if (b2 > nbins) b2 = ++b1;
  }
}

// [[Rcpp::export]]
DataFrame cts_to_mat(const DataFrame cts, int nbins, double dispersion, std::vector<double>& phi)
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
  
  
  cts_to_mat_core(N, &cts_bin1[0], &cts_bin2[0], &count[0], &lmu_nosig[0], &weight[0], nbins, dispersion, &phi[0],
                  &phihat[0], &phihat_var[0], &ncounts[0], &bin1[0], &bin2[0]);
  
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

RcppExport SEXP wgfl_perf(const DataFrame cts, double dispersion, int niter, int nbins,
        int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
        double lam,  double alpha, double inflate, int maxsteps, double converge)
{
    std::vector<int> trails_r = Rcpp::as<std::vector<int> >(trails_i);
    std::vector<int> breakpoints_r = Rcpp::as<std::vector<int> >(breakpoints_i);
    std::vector<double> z_r(breakpoints_r[ntrails-1], 0);
    std::vector<double> u_r(breakpoints_r[ntrails-1], 0);
    
    const int N = nbins*(nbins+1)/2; //size of fused lasso problem
    std::vector<double> phi(N, 0);
    printf("Fused lasso perf iteration with %d coefficients\n",phi.size());
    
    for (int step=0; step<niter; ++step) {
      const DataFrame mat = cts_to_mat(cts, nbins, dispersion, phi);
      std::vector<double> y_r = Rcpp::as<std::vector<double> >(mat["phihat"]);
      std::vector<double> w_r = Rcpp::as<std::vector<double> >(mat["weight"]);
      
      int res;
      res = graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails, &trails_r[0], &breakpoints_r[0],
                                         lam, alpha, inflate, maxsteps, converge,
                                         &phi[0], &z_r[0], &u_r[0]);
    }
    return wrap(phi);
}




RCPP_MODULE(gfl){
  using namespace Rcpp ;
 
  function("weighted_graphfl" , &weighted_graphfl  , "documentation for weighted_graphfl ");
  function("cts_to_mat" , &cts_to_mat  , "documentation for cts_to_mat ");
  //function("cts_to_mat_core" , &cts_to_mat_core  , "documentation for cts_to_mat_core ");
  function("wgfl_perf" , &wgfl_perf  , "documentation for wgfl_perf ");
} 

