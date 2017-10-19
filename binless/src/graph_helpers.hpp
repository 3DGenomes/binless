#ifndef CSNORM_GRAPH_TRAILS_HPP
#define CSNORM_GRAPH_TRAILS_HPP

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

struct Coordinate {
    int index,bin1,bin2;
};
typedef boost::adjacency_list
<boost::vecS, boost::vecS, boost::undirectedS, Coordinate> Graph;

Graph build_patch_graph(int nbins, double tol_val, const IntegerVector& bin1,
                        const IntegerVector& bin2, const NumericVector& value);

//mat must be sorted by bin1 and bin2 and will not be checked for that
IntegerVector get_patch_numbers(int nbins, double tol_val, const IntegerVector& bin1,
                                const IntegerVector& bin2, const NumericVector& value);

#endif
