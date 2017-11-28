#' @include binless.R
NULL

#' Fetch CSgroup indices from CSnorm object
#'
#' @param resolution 
#' @param group
#' @param raise boolean. If T raise an exception, otherwise return -1.
#' @param cs CSnorm object
#'
#' @return CSgroup object index
#' @keywords internal
#' @export
#'
#' @examples
get_cs_group_idx = function(cs, resolution, group, raise=T) {
  if(length(cs@groups)>0) {
    for (i in 1:length(cs@groups)) {
      if (cs@groups[[i]]@resolution==resolution & cs@groups[[i]]@group==group) return(i)
    }
  }
  if (raise==T) {
    stop("CSgroup not found. You must first call bin_all_datasets or group_datasets")
  } else {
    return(-1)
  }
}

#' Fetch CSinter object from CSgroup object
#'
#' @param csg
#' @param type 
#' @param threshold 
#' @param ref 
#' @param raise 
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
get_cs_interaction_idx = function(csg, type, threshold=-1, ref=NULL, raise=T) {
  if(length(csg@interactions)>0) {
    for (i in 1:length(csg@interactions)) {
      if (class(csg@interactions[[i]])!=type) next
      if (!is.null(ref)) if (csg@interactions[[i]]@ref!=ref) next
      if (threshold>=0) if (csg@interactions[[i]]@threshold!=threshold) next
      return(i)
    }
  }
  if (raise==T) {
    stop("CSinter not found. You must first detect the corresponding interactions and differences")
  } else {
    return(-1)
  }
}

#' Fetch grouped matrices from CSnorm object
#' 
#' 
#' @param cs CSnorm object.
#' @param resolution The resolution of the matrix, in bases
#' @param group character. Either "all" to return individual matrices, or any of 
#'   "condition", "replicate" or "enzyme" to return grouped matrices
#'   
#' @return a data.table containing the binned matrices
#' @export
#' 
#' @examples
get_binned_matrices = function(cs, resolution=cs@settings$base.res, group="all") {
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csb=cs@groups[[idx1]]
  return(csb@mat)
}

#' Fetch detected interactions from CSnorm object
#' 
#' @keywords internal
#' 
#' @examples
get_interactions = function(cs, type, resolution, group, threshold=-1, ref=NULL) {
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  idx2=get_cs_interaction_idx(csg, type, threshold, ref)
  mat=copy(csg@interactions[[idx2]]@mat)
  mat = merge(csg@mat,mat,by=c("name","bin1","bin2"))
  return(mat)
}

#' Fetch binned interactions from CSnorm object
#' 
#' @inheritParams get_binned_matrices
#' 
#' @return a data.table containing the interactions
#' @export
#' 
#' @examples
get_binned_interactions = function(cs, resolution=cs@settings$base.res, group="all", threshold=0.95) {
  type="CSsig"
  ref=NULL
  mat=get_interactions(cs, type, resolution, group, threshold, ref)
  return(mat)
}

#' Fetch binless interactions from CSnorm object
#' 
#' @inheritParams get_binned_matrices
#' 
#' @return a data.table containing the interactions
#' @export
#' 
#' @examples
get_binless_interactions = function(cs, resolution=cs@settings$base.res, group="all") {
  type="CSbsig"
  ref=NULL
  threshold=-1
  mat=get_interactions(cs, type, resolution, group, threshold, ref)
  mat[,c("signal","binless"):=list(exp(phi),exp(phi)*decaymat)]
  return(mat)
}

#' Fetch binned differences from CSnorm object
#' 
#' @inheritParams get_binned_matrices
#' 
#' @return a data.table containing the differences
#' @export
#' 
#' @examples
get_binned_differences = function(cs, ref, resolution=cs@settings$base.res, group="all", threshold=0.95) {
  type="CSdiff"
  mat=get_interactions(cs, type, resolution, group, threshold, ref)
  return(mat)
}

#' Fetch binless differences from CSnorm object
#' 
#' @inheritParams get_binned_matrices
#' 
#' @return a data.table containing the differences
#' @export
#' 
#' @examples
get_binless_differences = function(cs, ref, resolution=cs@settings$base.res, group="all") {
  type="CSbdiff"
  threshold=-1
  mat=get_interactions(cs, type, resolution, group, threshold, ref)
  mat[,difference:=exp(delta)]
  return(mat)
}

