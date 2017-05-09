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

// [[Rcpp::export]]
DataFrame cts_to_mat(const DataFrame cts, int nbins, double dispersion, const std::vector<double>& phi)
{
    //inputs
    size_t N = cts.nrows();
    IntegerVector cts_bin1 = cts["bin1"];
    IntegerVector cts_bin2 = cts["bin2"];
    NumericVector count = cts["count"];
    NumericVector lmu_nosig = cts["lmu.nosig"];
    NumericVector weight = cts["weight"];

    //outputs
    size_t nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
    NumericVector phihat(nbetas); //vectorized form
    NumericVector phihat_var(nbetas);
    NumericVector ncounts(nbetas);
    IntegerVector bin1(nbetas);
    IntegerVector bin2(nbetas);

    //walk through cts
    for (int i = 0; i < N; ++i) {
        //pos(b1,b2) starts at 0 and is the index in the 2D triangle grid
        //of the cell at coordinates (b1,b2), where b1 and b2 start at 1 and b2>=b1
        int b1 = cts_bin1(i);
        int b2 = cts_bin2(i);
        int pos = (b1-1)*(nbins+1) - (b1-1)*b1/2 + b2 - b1;
        double mu = std::exp( lmu_nosig(i) + phi[pos] );
        double z = count(i)/mu-1;
        double var = 1/mu + 1/dispersion;
        double w2v = weight(i)/(2*var);
        ncounts(pos) += weight(i);
        phihat_var(pos) += w2v;
        phihat(pos) += (z+phi[pos])*w2v;
    }

    //finish mat
    int b1 = 1;
    int b2 = 1;
    for (int i = 0; i < nbetas; ++i) {
        phihat_var(i) = 1/phihat_var(i);
        phihat(i) = phihat(i)*phihat_var(i);
        bin1(i) = b1;
        bin2(i) = b2++;
        if (b2 > nbins) b2 = ++b1;
    }

    bin1.attr("levels") = cts_bin1.attr("levels");
    bin2.attr("levels") = cts_bin2.attr("levels");
    bin1.attr("class") = CharacterVector::create("ordered", "factor");
    bin2.attr("class") = CharacterVector::create("ordered", "factor");

    return DataFrame::create(_["bin1"]=bin1, _["bin2"]=bin2, _["phihat"]=phihat,
            _["phihat.var"]=phihat_var, _["ncounts"]=ncounts, _["weight"]=1/phihat_var,
            _["diag.idx"]=bin2-bin1);
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
  function("wgfl_perf" , &wgfl_perf  , "documentation for wgfl_perf ");
} 

