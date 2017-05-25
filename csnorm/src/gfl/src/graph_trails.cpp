#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <assert.h>

typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

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

