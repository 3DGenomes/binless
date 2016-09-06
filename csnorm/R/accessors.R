#' @include csnorm.R
NULL

#' Fetch CSbinned indices from CSnorm object
#'
#' @param cs CSnorm object
#' @param resolution, dispersion.type see \code{\link{bin_all_datasets}}
#'
#' @return CSbinned object index
#' @keywords internal
#' @export
#'
#' @examples
get_cs_binned_idx = function(cs, resolution, dispersion.type) {
  stopifnot(length(cs@binned)>0)
  for (i in 1:length(cs@binned)) {
    if (cs@binned[[i]]@resolution==resolution && cs@binned[[i]]@dispersion.type==dispersion.type) return(i)
  }
}

#' Fetch binned matrix indices from CSnorm object
#'
#' @param cs CSnorm object
#' @param resolution,type,dispersion.type,dispersion.fun see 
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to 
#'   identify the input matrices.
#'
#' @return two integers
#' @keywords internal
#' @export
#'
#' @examples
get_matrices_idx = function(cs, resolution, type, dispersion.type, dispersion.fun) {
  stopifnot(length(cs@binned)>0)
  idx1=which(sapply(cs@binned, function(x){x@resolution==resolution && x@dispersion.type==dispersion.type}))
  stopifnot(length(idx1)==1)
  csb=cs@binned[[idx1]]
  if (type=="all") {
    idx2=0
  } else {
    idx2=which(sapply(csb@metadata[,.I],
                    function(x){csb@metadata[x,type]==type && csb@metadata[x,dispersion.fun]==deparse(dispersion.fun)}))-1
  }
  stopifnot(length(idx2)==1)
  return(c(idx1,idx2))
}

#' Fetch binned matrix from CSnorm object
#'
#' @inheritParams get_matrices_idx
#'
#' @return a data.table containing the binned matrices
#' @export
#'
#' @examples
get_matrices = function(cs, resolution, type, dispersion.type, dispersion.fun) {
  idx=get_matrices_idx(cs,resolution,type,dispersion.type,dispersion.fun)
  csb=cs@binned[[idx[1]]]
  if (idx[2]==0) return(csb@individual) else return(csb@grouped[[idx[2]]])
}

#' Fetch fit parameters from CSnorm object
#'
#' @param cs CSnorm object
#'
#' @return list of fit parameters
#' @export
#'
#' @examples
get_cs_parameters = function(cs) cs@par