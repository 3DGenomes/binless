#ifndef GRAPH_TRAILS_H
#define GRAPH_TRAILS_H

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>

typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;


std::vector<std::vector<int> > boost_triangle_grid_chain(int nrow);

List boost_chains_to_trails(const std::vector<std::vector<int> >& chains);

Graph build_2d_connectivity_graph(int nrow);

void print_2d_connectivity_graph(int nbins);

#endif

