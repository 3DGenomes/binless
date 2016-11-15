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
#' @param ref 
#' @param raise 
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
get_cs_interaction_idx = function(csm, type, ref, raise=T) {
  if(length(csm@interactions)>0) {
    for (i in 1:length(csm@interactions)) {
      if (csm@interactions[[i]]@type==type
          && csm@interactions[[i]]@ref==ref) return(i)
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
#' @return a data.table containing the interactions
#' @export
#' 
#' @examples
get_interactions = function(cs, type, resolution, group, ref) {
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  mat=csm@mat
  idx3=get_cs_interaction_idx(csm, type, ref)
  int=csm@interactions[[idx3]]@mat
  if (type=="interactions") {
    ret = merge(mat[,.(name,bin1,begin1,end1,bin2,begin2,end2,observed,signal,signal.sd,normalized)],
                int[,.(name,bin1,bin2,signal.signif=signal,direction)], by=c("name","bin1","bin2"))
    return(ret)
  } else {
    ret=merge(ret, mat[name==ref,.(bin1,bin2,ref.observed=observed,ref.expected=expected,
                                   ref.normalized=normalized,ref.signal=signal)], by=c("bin1","bin2"))
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
  foreach (gen=cs@design[,unique(genomic)], .combine=rbind) %do% {
    dsets=cs@design[genomic==gen,name]
    csnorm:::generate_genomic_biases(biases=cs@biases[name%in%dsets], beta_iota=cs@par$beta_iota[gen,],
                                     beta_rho=cs@par$beta_rho[gen,], bf_per_kb=cs@settings$bf_per_kb,
                                     points_per_kb = points_per_kb)[,.(pos,log_iota,log_rho,genomic=gen)]
  }
}

#' update diagnostics data table or create it if not existing
#'
#' @return
#' @keywords internal
#'
#' @examples
update_diagnostics = function(cs, step, leg, out, runtime) {
  params=data.table(step=step,leg=leg,out=list(out),out.last=tail(out,n=1),value=cs@par$value,runtime=runtime)
  tmp=as.data.table(lapply(cs@par,list))
  #remove entries that are too heavy and redundant
  if ("log_decay" %in% names(tmp)) tmp[,log_decay:=NULL]
  if ("biases" %in% names(tmp)) tmp[,biases:=NULL]
  #merge with previous
  params=cbind(params,tmp)
  if (is.data.table(cs@diagnostics$params)) params=rbind(cs@diagnostics$params,params,fill=T)
  return(params)
}


#' Returns one param's values during optimization
#'
#' @keywords internal
#' @export
#'
#' @examples
get_all_values = function(cs, param, trans) {
  #get value in tmp as vector of lists, remove NULL lists
  values=cs@diagnostics$params[,.(step=step+((.I-1)%%3)/3,leg,tmp=get(param))][!sapply(tmp,is.null)]
  #melt it
  melted=as.data.table(values[,melt(tmp)])
  if ("Var1" %in% names(melted)) {
    if ("Var2" %in% names(melted)) {
      melted[,variable:=paste0(param,".",Var1,".",Var2)]
    } else {
      melted[,variable:=paste0(param,".",Var1)]
    }
  } else {
    melted[,variable:=param]
  }
  #merge it back
  values[,L1:=.I]
  melted=merge(values,melted,by="L1")[,.(variable,step,leg,value)]
  #transform if necessary
  if (!is.na(trans)) {
    if (trans=="log") {
      melted[,c("variable","value"):=list(paste("log",variable),log(value))]
    } else if (trans=="exp") {
      melted[,c("variable","value"):=list(paste("exp",variable),exp(value))]
    } else {
      stop("unsupported transformation ",trans)
    }
  }
  #return
  setkey(melted,variable,step,leg)
  melted
}
