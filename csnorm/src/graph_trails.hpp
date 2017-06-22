#ifndef CSNORM_GRAPH_TRAILS_HPP
#define CSNORM_GRAPH_TRAILS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

struct Coordinate {
    int index,bin1,bin2;
};
typedef boost::adjacency_list
<boost::vecS, boost::vecS, boost::undirectedS, Coordinate> Graph;

std::vector<std::vector<int> > boost_triangle_grid_chain(int nrow);

List boost_chains_to_trails(const std::vector<std::vector<int> >& chains);

Graph build_patch_graph(int nrow, const DataFrame mat, double tol_val);

//mat must be sorted by bin1 and bin2 and will not be checked for that
List boost_build_patch_graph_components(int nbins, const DataFrame mat,
                                        double tol_val);

#endif
