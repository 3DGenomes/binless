#ifndef CSNORM_GRAPH_TRAILS_HPP
#define CSNORM_GRAPH_TRAILS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>

struct Coordinate { int index,bin1,bin2; };
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, Coordinate> Graph;

struct edge_within_patch {
  
  edge_within_patch() { }
  
  edge_within_patch(const Graph& G_, double* value_, //ptr to begin of std::vector<double> of size |E|
                    double tol_val_)
    : G(G_), value(value_), tol_val(tol_val_) { }
  
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return std::abs(get(value, boost::source(e, G))-get(value, boost::target(e, G))) < tol_val;
  }
  
  Graph G;
  double* value;
  double tol_val;
};

std::vector<std::vector<int> > boost_triangle_grid_chain(int nrow);

List boost_chains_to_trails(const std::vector<std::vector<int> >& chains);

Graph build_2d_connectivity_graph(int nrow);

std::vector<double> report_values_in_graph(Graph G, const DataFrame mat);

//mat must be sorted by bin1 and bin2 and will not be checked for that
List boost_build_patch_graph_components(int nbins, const DataFrame mat, double tol_val);

#endif
