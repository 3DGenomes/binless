#' Fast binless signal calculation
#' 
#' This function computes a fast approximation to the binless normalization of a given dataset
#' Most notably, it does not produce a statistically significant output, since the fusion and
#' significance thresholds must be set by the user. Other differences include a pointwise estimate
#' of the decay and genomic biases, one decay and one genomic bias for all datasets, and
#' a Poisson model of the counts (infinite dispersion). This implies that empty rows and
#' counter diagonals are not admissible and will throw an error. This fast model is meant for a quick
#' overview of the data, and should not be used for high-resolution and/or statistical analyses.
#'
#' @param obs DataFrame containing at least 4 named columns: name, bin1, bin2 and observed.
#'  The first three must be convertible to integer vectors (like factors) and must range from
#'  1 to N contiguously. For example, if there are 2 datasets, name should range from 1 to 2,
#'  and similarly for the bin columns bin1 and bin2. observed is the number of reads that fall
#'  at that coordinate. The provided matrices must be dense, i.e. report zero counts.
#' @param nbins unsigned integer referring to the number of bins in the provided matrix.
#' @param lam2 numeric positive value for the fusion penalty.
#' @param nouter unsigned The maximum number of iterations that should be performed (default 20)
#' @param tol_val double tolerance on the values for convergence and the fused lasso (default 1e-1)
#' @param bg_steps unsigned the maximum number of initial steps where no signal is fitted
#' 
#' @export
#' @name fast_binless 
NULL

#' Fast binless difference calculation
#' 
#' Once normalized with \code{\link{fast_binless}}, differences with respect to a reference can be computed.
#'
#' @ref unsigned integer corresponding to the index (starting at 1) of the dataset to use as reference
#' @inheritParams fast_binless
#' 
#' @export
#' @name fast_binless_difference 
NULL

#' Build binned matrix using raw cts data
#' 
#' @param cs 
#' @param cts 
#' @param resolution 
bin_cts_to_mat = function(cs, cts, resolution) {
  nbins = cts[,nlevels(bin1)]
  alpha = cs@par$alpha
  mat=foreach(n=cts[,unique(name)],.combine=rbind) %do% {
    metadata = binless:::get_signal_metadata(cs, cts[name==n], resolution)
    mat = as.data.table(binless:::rcpp_cts_to_signal_mat(nbins, alpha, cts[name==n], metadata))
    mat[,name:=n]
    mat
  }
  return(mat)
}




