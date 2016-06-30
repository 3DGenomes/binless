#' @include csnorm.R
NULL

#' Fetch CSbinned matrices from CSnorm object
#'
#' @param cs CSnorm object
#' @param idx index of the binned matrix
#'
#' @return CSbinned matrix
#' @export
#'
#' @examples
get_cs_binned = function(cs, idx, mat=c("CS","ICE","RAW")) {
  match.arg(mat)
  stopifnot(length(mat)==1)
  if (mat=="CS") {
    cs@binned[[idx]]@mat
  } else if (mat=="RAW") {
    cs@binned[[idx]]@raw
  } else {
    cs@binned[[idx]]@ice
  }
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