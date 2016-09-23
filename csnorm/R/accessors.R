#' @include csnorm.R
NULL

#' Fetch CSbinned indices from CSnorm object
#'
#' @param resolution 
#' @param raise boolean. If T raise an exception, otherwise return -1.
#' @param cs CSnorm object
#'
#' @return CSbinned object index
#' @keywords internal
#' @export
#'
#' @examples
get_cs_binned_idx = function(cs, resolution, raise=T) {
  if(length(cs@binned)>0) {
    for (i in 1:length(cs@binned)) {
      if (cs@binned[[i]]@resolution==resolution) return(i)
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
#' @param group 
#' @param raise boolean. If T raise an exception, otherwise return -1.
#'
#' @return CSmatrix object index
#' @keywords internal
#' @export
#'
#' @examples
get_cs_matrix_idx = function(csb, group, raise=T) {
  for (i in 1:length(csb@grouped)) {
    if (all(csb@grouped[[i]]@group==group)) return(i)
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
#' @param ref 
#' @param raise 
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
get_cs_interaction_idx = function(csm, type, threshold, ref, raise=T) {
  if(length(csm@interactions)>0) {
    for (i in 1:length(csm@interactions)) {
      if (csm@interactions[[i]]@type==type
          && csm@interactions[[i]]@threshold==threshold && csm@interactions[[i]]@ref==ref) return(i)
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
#' @param group character. Either "all" to return individual matrices, or any of 
#'   "condition", "replicate" or "enzyme" to return grouped matrices
#'   
#' @return a data.table containing the binned matrices
#' @export
#' 
#' @examples
get_matrices = function(cs, resolution, group) {
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  return(csb@grouped[[idx2]]@mat)
}


#' Fetch detected interactions from CSnorm object
#' 
#' @param type character. Either "interactions" or "differences"
#' @inheritParams get_matrices
#' @param ref character. The 
#' @param threshold 
#' @return a data.table containing the interactions
#' @export
#' 
#' @examples
get_interactions = function(cs, type, resolution, group,
                            ref, threshold) {
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  mat=csm@mat
  idx3=get_cs_interaction_idx(csm, type, threshold, ref)
  int=csm@interactions[[idx3]]@mat
  ret=merge(mat,int,by=c("name","bin1","bin2"))
  if (type=="interactions") {
    return(ret)
  } else {
    ret=merge(ret, mat[name==ref,.(bin1,bin2,ref.observed=observed,ref.expected=expected,
                                   ref.normalized=normalized,ref.icelike=icelike)], by=c("bin1","bin2"))
    return(ret)
  }
}

#' Predict model parameters on a subset of the input data
#'
#' @param cs 
#' @param ncounts 
#'
#' @return
#' @export
#'
#' @examples
get_predicted_subset = function(cs, ncounts=100000) {
  if (length(cs@par)==0) stop("You should first normalize the datasets")
  counts=cs@counts[sample(.N,min(.N,ncounts))]
  setkey(counts,name,id1,id2)
  mat=csnorm_predict_all(cs, counts, verbose=F)
  return(mat)
}

#' Predict genomic biases at evenly spaced positions across the genome
#'
#' @param cs 
#' @param points_per_kb 
#'
#' @return
#' @export
#'
#' @examples
get_genomic_biases = function(cs, points_per_kb=10) {
  if (length(cs@par)==0) stop("You should first normalize the datasets")
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,dset=j)]
}
