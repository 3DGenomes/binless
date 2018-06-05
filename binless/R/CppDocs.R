#' Fast binless signal calculation
#' 
#' This function computes a fast approximation to the binless normalization of a given dataset
#' Most notably, it does not produce a statistically significant output, since the fusion and
#' significance thresholds must be set by the user. Other differences include a pointwise estimate
#' of the decay and genomic biases, and one decay and one genomic bias for all datasets. This fast 
#' model is meant for a quick overview of the data.
#'
#' @usage fast_binless(obs, nbins, alpha, lam2, lam1 = 0, nouter = 25, tol_val = 2e-1,
#'  bg_steps = 5, free_decay = 10000)
#'
#' @param obs DataFrame containing at least 4 named columns: name, bin1, bin2 and observed.
#'  The first three must be convertible to integer vectors (like factors) and must range from
#'  1 to N contiguously. For example, if there are 2 datasets, name should range from 1 to 2,
#'  and similarly for the bin columns bin1 and bin2. observed is the number of reads that fall
#'  at that coordinate. The provided matrices must be dense, i.e. report zero counts.
#' @param nbins unsigned integer referring to the number of bins in the provided matrix.
#' @param alpha the dispersion of the negative binomial
#' @param lam2 numeric positive value for the fusion penalty, or a vector of
#' the same size as the number of datasets.
#' @param lam1 numeric positive value (or vector) for the significance threshold (default is zero)
#' @param nouter unsigned The maximum number of iterations that should be performed (default 25)
#' @param tol_val double tolerance on the values for convergence and the fused lasso (default 2e-1)
#' @param bg_steps unsigned the maximum number of initial steps where no signal is fitted (default 5)
#' @param free_decay the distance in bases up to which the decay is not forced to decrease (default 10000)
#' 
#' @export
#' @name fast_binless 
NULL

#' Fast binless difference calculation
#' 
#' Once normalized with \code{\link{fast_binless}}, differences with respect to a reference can be computed.
#' 
#' @usage fast_binless_difference(out, ref, alpha, lam2, lam1 = 0, tol_val = 2e-1)
#'
#' @param ref unsigned integer corresponding to the index (starting at 1) of the dataset to use as reference
#' @param lam2 numeric positive value for the fusion penalty, or a vector of
#' the same size as the number of datasets minus one.
#' @param lam1 numeric positive value (or vector of size n-1) for the significance threshold (default is zero)
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


#' Generate an empty matrix in triangular form.
#' 
#' arguments are taken as ordered factors
#' returned matrix is of size n_names * n_bins * (n_bins+1) / 2
#' n_names and n_bins are the number of levels in the name and bins factors
#'
#' @param name an ordered factor of names
#' @param bins an ordered factor of bins
#'
#' @return
#'
#' @examples
create_empty_matrix = function(name, bins) {
  dt=binless:::create_empty_matrix_cpp(name,bins)
  dt=setDT(dt,key=c("name","bin1","bin2"))
  return(dt)
}

