#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>

#include "graph_trails.hpp"

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

Graph build_patch_graph(int nrow, const DataFrame mat, double tol_val) {
  //create graph
  int ntotal = nrow*(nrow+1)/2;
  Graph G(ntotal);
  //std::cout << "Graph with " << ntotal << " vertices" << std::endl;
  //set coordinates in the 2d plane
  Graph::vertex_descriptor v = *boost::vertices(G).first;
  for (int b1=1, i=0; b1<=nrow; ++b1) {
    for (int b2=b1; b2<=nrow; ++b2, ++v, ++i) {
      G[v].index = i;
      G[v].bin1 = b1;
      G[v].bin2 = b2;
      //std::cout << "vertex " << v << " has coordinates (" << b1 << "," << b2 << ")" << std::endl;
    }
  }
  //retrieve data from matrix
  IntegerVector bin1 = mat["bin1"];
  IntegerVector bin2 = mat["bin2"];
  NumericVector value = mat["value"];
  //add edges in 2d triangle grid if they are in the same patch
  v = *boost::vertices(G).first+1; //start at 2nd vertex
  for (int l=nrow-1; l>=1; ++v, --l) {
    for (int j=1; j<=l; ++v, ++j) {
      int i=G[v].index;
      if (G[v].bin1 != bin1(i) || G[v].bin2 != bin2(i)) {
        throw std::invalid_argument("mat must be ordered!");
      }
      if (std::abs(value(i-1)-value(i)) < tol_val) boost::add_edge(v-1,v,G);
      if (std::abs(value(i)-value(i+l)) < tol_val) boost::add_edge(v,v+l,G);
    }
  }
  return(G);
}


//mat must be sorted by bin1 and bin2 and will not be checked for that
List boost_build_patch_graph_components(int nbins, const DataFrame mat, double tol_val) {
  //build graph with edges only between vertices with equal values
  Graph fG = build_patch_graph(nbins, mat, tol_val);
  //deduce connected components
  std::vector<int> component(boost::num_vertices(fG));
  int num = boost::connected_components(fG, &component[0]);
  /*std::cout << "Total number of components: " << num << std::endl;
  std::vector<int>::size_type i;
  for (i = 0; i != component.size(); ++i)
    std::cout << "Vertex " << i <<" is in component " << component[i] << std::endl;
  std::cout << std::endl;*/
  return List::create(_["no"]=num, _["membership"]=component);
}

