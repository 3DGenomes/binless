#' @include csnorm.R
NULL

#' Fetch CSbinned indices from CSnorm object
#'
#' @param resolution 
#' @param dispersion.type 
#' @param raise boolean. If T raise an exception, otherwise return -1.
#' @param cs CSnorm object
#'
#' @return CSbinned object index
#' @keywords internal
#' @export
#'
#' @examples
get_cs_binned_idx = function(cs, resolution, dispersion.type, raise=T) {
  if(length(cs@binned)>0) {
    for (i in 1:length(cs@binned)) {
      if (cs@binned[[i]]@resolution==resolution && cs@binned[[i]]@dispersion.type==dispersion.type) return(i)
    }
  }
  if (raise==T) {
    stop("CSbinned not found. You must first call bin_all_datasets")
  } else {
    return(-1)
  }
}

#' Fetch CSmatrix indices from CSbinned object
#'
#' @param csb CSbinned object
#' @param type 
#' @param dispersion.fun 
#' @param raise boolean. If T raise an exception, otherwise return -1.
#'
#' @return CSmatrix object index
#' @keywords internal
#' @export
#'
#' @examples
get_cs_matrix_idx = function(csb, type, dispersion.fun, raise=T) {
  for (i in 1:length(csb@grouped)) {
    if (csb@grouped[[i]]@type==type && csb@grouped[[i]]@dispersion.fun==deparse(dispersion.fun,nlines=1)) return(i)
  }
  if (raise==T) {
    stop("CSmatrix not found. You must first call group_datasets")
  } else {
    return(-1)
  }
}

#' Fetch CSinter object from CSmatrix object
#'
#' @param csm 
#' @param type 
#' @param threshold 
#' @param normal.approx 
#' @param ref 
#' @param raise 
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
get_cs_interaction_idx = function(csm, type, threshold, normal.approx, ref, raise=T) {
  if(length(csm@interactions)>0) {
    for (i in 1:length(csm@interactions)) {
      if (csm@interactions[[i]]@type==type && csm@interactions[[i]]@threshold==threshold
          && csm@interactions[[i]]@normal.approx==normal.approx && csm@interactions[[i]]@ref==ref) return(i)
    }
  }
  if (raise==T) {
    stop("CSinter not found. You must first call detect_interactions or detect_differences")
  } else {
    return(-1)
  }
}

#' Fetch binned matrix from CSnorm object
#' 
#' 
#' @param cs CSnorm object.
#' @param resolution The resolution of the matrix, in bases
#' @param type character. Either "all" to return individual matrices, or any of 
#'   "condition", "replicate" or "enzyme" to return grouped matrices
#' @param dispersion.type,dispersion.fun dispersion parameters, as provided to
#'   the binning and grouping funtions.
#'   
#' @return a data.table containing the binned matrices
#' @export
#' 
#' @examples
get_matrices = function(cs, resolution, type, dispersion.type, dispersion.fun) {
  idx1=get_cs_binned_idx(cs, resolution, dispersion.type, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, type, dispersion.fun, raise=T)
  return(csb@grouped[[idx2]]@mat)
}


#' Fetch detected interactions from CSnorm object
#' 
#' @param interaction.type character. Either "interactions" or "differences"
#' @inheritParams get_matrices
#' @param ref 
#' @param threshold 
#' @param normal.approx 
#' @return a data.table containing the interactions
#' @export
#' 
#' @examples
get_interactions = function(cs, interaction.type, resolution, type, dispersion.type,
                            dispersion.fun, ref, threshold, normal.approx) {
  idx1=get_cs_binned_idx(cs, resolution, dispersion.type, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, type, dispersion.fun, raise=T)
  csm=csb@grouped[[idx2]]
  idx3=get_cs_interaction_idx(csm, interaction.type, threshold, normal.approx, ref)
  return(csm@interactions[[idx3]]@mat)
}
