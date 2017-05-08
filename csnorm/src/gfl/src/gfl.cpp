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
NumericVector phi_to_cts(const DataFrame cts, int nbins, const std::vector<double>& phi)
{
  NumericMatrix phimat(nbins, nbins);
  int pos = 0;
  for (int i = 0; i < nbins; ++i) {
    for (int j = i; j < nbins; ++j) {
      phimat(i,j) = phi[pos];
      pos++;
      if (pos == phi.size()) break;
    }
    if (pos == phi.size()) break;
  }
  
  IntegerVector bin1 = cts["bin1"];
  IntegerVector bin2 = cts["bin2"];
  NumericVector phivec(cts.nrows());
  for(int i=0; i<cts.nrows(); ++i)
    phivec(i) = phimat(bin1[i]-1, bin2[i]-1);
  
  return(phivec);
}

// [[Rcpp::export]]
DataFrame cts_to_mat(const DataFrame cts, int nbins, double dispersion, const std::vector<double>& phi_mat)
{
    //inputs
    IntegerVector bin1 = cts["bin1"];
    IntegerVector bin2 = cts["bin2"];
    NumericVector count = cts["count"];
    NumericVector lmu_nosig = cts["lmu.nosig"];
    NumericVector phi = phi_to_cts(cts,nbins,phi_mat);
    NumericVector weight = cts["weight"];

    //outputs
    size_t N = cts.nrows();
    NumericMatrix phihat(nbins, nbins);
    NumericMatrix phihat_var(nbins, nbins);
    NumericMatrix ncounts(nbins, nbins);
    NumericVector phihat_r; //vectorized form
    NumericVector phihat_var_r;
    NumericVector ncounts_r;
    IntegerVector bin1_r;
    IntegerVector bin2_r;

    //walk through cts
    for (int i = 0; i < N; ++i) {
        double mu = std::exp( lmu_nosig(i) + phi(i) );
        double z = count(i)/mu-1;
        double var = 1/mu + 1/dispersion;
        double w2v = weight(i)/(2*var);
        ncounts(bin1(i)-1, bin2(i)-1) += weight(i);
        phihat_var(bin1(i)-1, bin2(i)-1) += w2v;
        phihat(bin1(i)-1, bin2(i)-1) += (z+phi(i))*w2v;
    }

    //finish mat
    for (int i = 0; i < nbins; ++i) {
        for (int j = i; j < nbins; ++j) {
            phihat_var(i,j) = 1/phihat_var(i,j);
            phihat(i,j) = phihat(i,j)*phihat_var(i,j);
        }
    }

    //vectorize and put into data frame
    for (int i = 0; i < nbins; ++i) {
        for (int j = i; j < nbins; ++j) {
            bin1_r.push_back(i+1);
            bin2_r.push_back(j+1);
            phihat_r.push_back(phihat(i,j));
            phihat_var_r.push_back(phihat_var(i,j));
            ncounts_r.push_back(ncounts(i,j));
        }
    }
    bin1_r.attr("levels") = bin1.attr("levels");
    bin2_r.attr("levels") = bin2.attr("levels");
    bin1_r.attr("class") = CharacterVector::create("ordered", "factor");
    bin2_r.attr("class") = CharacterVector::create("ordered", "factor");

    return DataFrame::create(_["bin1"]=bin1_r, _["bin2"]=bin2_r, _["phihat"]=phihat_r,
            _["phihat.var"]=phihat_var_r, _["ncounts"]=ncounts_r, _["weight"]=1/phihat_var_r,
            _["diag.idx"]=bin2_r-bin1_r);
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
  //function("cts_to_mat" , &cts_to_mat  , "documentation for cts_to_mat ");
  //function("phi_to_cts" , &phi_to_cts  , "documentation for phi_to_cts ");
  function("wgfl_perf" , &wgfl_perf  , "documentation for wgfl_perf ");
} 

