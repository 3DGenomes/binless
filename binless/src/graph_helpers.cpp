#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>

#include "graph_helpers.hpp"

Graph build_patch_graph(int nbins, double tol_val, const IntegerVector& bin1,
                        const IntegerVector& bin2, const NumericVector& value) {
    //create graph
    int ntotal = nbins*(nbins+1)/2;
    Graph G(ntotal);
    //std::cout << "Graph with " << ntotal << " vertices" << std::endl;
    //set coordinates in the 2d plane
    Graph::vertex_descriptor v = *boost::vertices(G).first;
    for (int b1=1, i=0; b1<=nbins; ++b1) {
        for (int b2=b1; b2<=nbins; ++b2, ++v, ++i) {
            G[v].index = i;
            G[v].bin1 = b1;
            G[v].bin2 = b2;
            //std::cout << "vertex " << v << " has coordinates (" << b1 << "," << b2 << ")" << std::endl;
        }
    }
    //add edges in 2d triangle grid if they are in the same patch
    v = *boost::vertices(G).first+1; //start at 2nd vertex
    for (int l=nbins-1; l>=1; ++v, --l) {
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
IntegerVector get_patch_numbers(int nbins, double tol_val, const IntegerVector& bin1,
                                const IntegerVector& bin2, const NumericVector& value) {
    //build graph with edges only between vertices with equal values
    Graph fG = build_patch_graph(nbins, tol_val, bin1, bin2, value); //assume the lasso has been computed at precision greater than tol_val
    //deduce connected components
    std::vector<int> component(boost::num_vertices(fG));
    int num = boost::connected_components(fG, &component[0]);
    /*std::cout << "Total number of components: " << num << std::endl;
    std::vector<int>::size_type i;
    for (i = 0; i != component.size(); ++i)
      std::cout << "Vertex " << i <<" is in component " << component[i] << std::endl;
    std::cout << std::endl;*/
    /*for (unsigned i=0; i<component.size(); ++i) Rcpp::Rcout << " i= " << i << " bin1= " << bin1(i) << " bin2= " << bin2(i) << " value= " << value(i) << " patchno= " << component[i] << "\n";*/
    return Rcpp::wrap(component);
}







