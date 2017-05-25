#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <assert.h>

#include "gfl_c.c"
#include "graph_fl.c"
#include "tf_dp.c"
#include "utils.c"



typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;


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

// [[Rcpp::export]]
DataFrame cts_to_diff_mat(const DataFrame cts, const DataFrame ref, int nbins, double dispersion,
                          std::vector<double>& phi_ref, std::vector<double>& delta, int diag_rm)
{
  const DataFrame mat_ref = cts_to_signal_mat(ref, nbins, dispersion, phi_ref, diag_rm);
  
  std::vector<double> phi_oth;
  phi_oth.reserve(delta.size());
  for (int i=0; i<delta.size(); ++i) phi_oth[i] = phi_ref[i] + delta[i];
  const DataFrame mat_oth = cts_to_signal_mat(cts, nbins, dispersion, phi_oth, diag_rm);
  
  IntegerVector bin1 = mat_ref["bin1"];
  IntegerVector bin2 = mat_ref["bin2"];
  NumericVector phihat_ref = mat_ref["phihat"];
  NumericVector phihat_var_ref = mat_ref["phihat.var"];
  NumericVector phihat = mat_oth["phihat"];
  NumericVector phihat_var = mat_oth["phihat.var"];
  NumericVector ncounts_oth = mat_oth["ncounts"];
  NumericVector ncounts_ref = mat_ref["ncounts"];
  NumericVector deltahat = phihat-phihat_ref;
  NumericVector deltahat_var = phihat_var+phihat_var_ref;
  NumericVector ncounts = ncounts_ref+ncounts_oth;
  NumericVector weight = 1/deltahat_var;
  NumericVector didx = mat_ref["diag.idx"];
  
  return DataFrame::create(_["bin1"]=bin1, _["bin2"]=bin2,
                           _["phihat"]=phihat, _["phihat.var"]=phihat_var,
                           _["phihat.ref"]=phihat_ref, _["phihat.var.ref"]=phihat_var_ref,
                           _["deltahat"]=deltahat, _["deltahat.var"]=deltahat_var, _["ncounts"]=ncounts,
                             _["weight"]=weight, _["diag.idx"]=didx);
}

// [[Rcpp::export]]
List wgfl_diff_perf_warm(const DataFrame cts, const DataFrame ref, double dispersion, int nouter, int nbins,
                         int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                         double lam,  double alpha, double inflate, int ninner, double converge,
                         int diag_rm, NumericVector z_i, NumericVector u_i, NumericVector phi_ref_i,
                         NumericVector delta_i)
{
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  std::vector<int> trails_r = as<std::vector<int> >(trails_i);
  std::vector<int> breakpoints_r = as<std::vector<int> >(breakpoints_i);
  std::vector<double> z_r = as<std::vector<double> >(z_i);
  std::vector<double> u_r = as<std::vector<double> >(u_i);
  std::vector<double> phi_ref_r = as<std::vector<double> >(phi_ref_i);
  std::vector<double> delta_r = as<std::vector<double> >(delta_i);
  std::vector<double> delta_old = delta_r;
  
  int step;
  int res=0;
  //printf(" Perf iteration: start with lam=%f alpha=%f phi_ref[0]=%f delta[0]=%f z[0]=%f u[0]=%f\n",
  //       lam, alpha, phi_ref_r[0], delta_r[0], z_r[0], u_r[0]);
  for (step=1; step<=nouter; ++step) {
    const DataFrame mat = cts_to_diff_mat(cts, ref, nbins, dispersion, phi_ref_r, delta_r, diag_rm);
    std::vector<double> phihat_ref = Rcpp::as<std::vector<double> >(mat["phihat.ref"]);
    std::vector<double> phihat_var_ref = Rcpp::as<std::vector<double> >(mat["phihat.var.ref"]);
    std::vector<double> phihat = Rcpp::as<std::vector<double> >(mat["phihat"]);
    std::vector<double> phihat_var = Rcpp::as<std::vector<double> >(mat["phihat.var"]);
    std::vector<double> y_r = Rcpp::as<std::vector<double> >(mat["deltahat"]);
    std::vector<double> w_r = Rcpp::as<std::vector<double> >(mat["weight"]);
    
    res += graph_fused_lasso_weight_warm (N, &y_r[0], &w_r[0], ntrails, &trails_r[0], &breakpoints_r[0],
                                          lam, &alpha, inflate, ninner, converge,
                                          &delta_r[0], &z_r[0], &u_r[0]);
    for (int i=0; i<N; ++i) {
      if (phihat_var_ref[i]==INFINITY && phihat_var[i]==INFINITY) {
        phi_ref_r[i]=(phihat_ref[i]+phihat[i])/2;
      } else {
        phi_ref_r[i] = (phihat_ref[i]/phihat_var_ref[i] + (phihat[i]-delta_r[i])/phihat_var[i])
        /(1/phihat_var_ref[i]+1/phihat_var[i]);
      }
    }
    
    double maxval = std::abs(delta_r[0]-delta_old[0]);
    for (int i=1; i<N; ++i) maxval = std::max(std::abs(delta_r[i]-delta_old[i]), maxval);
    //printf(" Iteration %d with alpha=%f reached maxval=%.5e after %d steps (phi_ref[0]=%f delta[0]=%f z[0]=%f u[0]=%f)\n",
    //       step,alpha,maxval,res,phi_ref_r[0],delta_r[0], z_r[0], u_r[0]);
    if (maxval<converge) break;
    delta_old = delta_r;
  }
  //printf(" Perf iteration: end   with lam=%f alpha=%f phi_ref[0]=%f delta[0]=%f z[0]=%f u[0]=%f nouter=%d ninner=%d\n",
  //       lam, alpha, phi_ref_r[0], delta_r[0], z_r[0], u_r[0], step, res);
  return List::create(_["phi.ref"]=wrap(phi_ref_r), _["delta"]=wrap(delta_r), _["alpha"]=wrap(alpha),
                      _["mat"]=cts_to_diff_mat(cts, ref, nbins, dispersion, phi_ref_r, delta_r, diag_rm),
                      _["z"]=wrap(z_r), _["u"]=wrap(u_r), _["nouter"]=step, _["ninner"]=res);
}

// [[Rcpp::export]]
List wgfl_diff_perf(const DataFrame cts, const DataFrame ref, double dispersion, int nouter, int nbins,
                    int ntrails, const NumericVector trails_i, const NumericVector breakpoints_i,
                    double lam,  double alpha, double inflate, int ninner, double converge, int diag_rm)
{
  NumericVector z_i(breakpoints_i(ntrails-1));
  NumericVector u_i(breakpoints_i(ntrails-1));
  const int N = nbins*(nbins+1)/2; //size of fused lasso problem
  NumericVector phi_ref_i(N);
  NumericVector delta_i(N);
  //printf("Fused lasso cold perf iteration with %d coefficients\n",phi.size());
  return wgfl_diff_perf_warm(cts, ref, dispersion, nouter, nbins, ntrails, trails_i, breakpoints_i,
                             lam, alpha, inflate, ninner, converge, diag_rm, z_i, u_i, phi_ref_i, delta_i);
}

// [[Rcpp::export]]
std::vector<std::vector<int> > boost_triangle_grid_chain(int nrow) {
  int ntotal = nrow*(nrow+1)/2-1;
  std::vector<std::vector<int> > chains;
  int l = nrow;
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
  for (int U1=2; U1<=nrow; ++U1) {
    int Ui=U1;
    current = std::vector<int>(1,Ui-1);
    for (int i=1; i<U1; ++i) {
      int Uip1 = Ui + nrow - i;
      current.push_back(Uip1-1);
      Ui=Uip1;
    }
    chains.push_back(current);
  }
  return(chains);
}

// [[Rcpp::export]]
List boost_chains_to_trails(const std::vector<std::vector<int> >& chains) {
  std::vector<int> trails;
  std::vector<int> breakpoints;
  for (std::vector<std::vector<int> >::const_iterator it = chains.begin() ; it != chains.end(); ++it) {
    if (trails.size()>0) breakpoints.push_back(trails.size());
    trails.insert(trails.end(), it->begin(), it->end());
  }
  if (trails.size()>0) breakpoints.push_back(trails.size());
  return(List::create(_["ntrails"]=breakpoints.size(), _["trails"]=trails, _["breakpoints"]=breakpoints));
}

Graph build_2d_connectivity_graph(int nrow) {
  int ntotal = nrow*(nrow+1)/2-1;
  Graph G(ntotal+1);
  for (int i=1, l=nrow-1; l>=1; ++i, --l) {
    for (int j=1; j<=l; ++i, ++j) {
      boost::add_edge(i-1,i,G);
      boost::add_edge(i,i+l,G);
      //std::cout << "added edge " << i << " - " << i-1 << std::endl;
      //std::cout << "added edge " << i << " - " << i+l << std::endl;
    }
  }
  return(G);
}

// [[Rcpp::export]]
void print_2d_connectivity_graph(int nbins) {
  
  Graph G = build_2d_connectivity_graph(nbins);
  
  std::vector<int> component(num_vertices(G));
  int num = boost::connected_components(G, &component[0]);
  
  std::vector<int>::size_type i;
  std::cout << "Total number of components: " << num << std::endl;
  for (i = 0; i != component.size(); ++i)
    std::cout << "Vertex " << i <<" is in component " << component[i] << std::endl;
  std::cout << std::endl;
}

