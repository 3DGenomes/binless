#' @include binless.R
NULL

#' Perform binning of data
#' 
#' fast IRLS and zero counts approximation
#'
#' @param cts cts
#' @param dispersion dispersion
#' @param verbose boolean.
#'
#' @return " "
#' @keywords internal
#' @export
#'
#' @examples " "
predict_binned_matrices_irls = function(cts, dispersion) {
  #predict means
  cts[,c("decay","biases","mu.nosig"):=list(exp(log_decay),exp(log_bias),exp(lmu.nosig))]
  mat=cts[,.(observed=sum(count*nobs),
             nobs=sum(nobs),
             biasmat=weighted.mean(biases,nobs),
             decaymat=weighted.mean(decay,nobs),
             background=sum(mu.nosig*nobs),
             background.sd=sqrt(sum((mu.nosig+mu.nosig^2/dispersion)*nobs)),
             residual=sum(count*nobs)/sum(mu*nobs),
             normalized=sum(count*nobs)/sum(exp(lmu.nosig-log_decay)*nobs))
          ,keyby=c("name","bin1","bin2")]
  #fill missing data
  bl1=cts[,levels(bin1)]
  bins=ordered(1:length(bl1))
  levels(bins) <- bl1
  fullmat=create_empty_matrix(name=cts[,ordered(unique(name))], bins=bins)
  mat=mat[fullmat]
  mat[is.na(observed),c("nobs","observed","normalized"):=list(0,0,0)]
  mat[,diag.idx:=unclass(bin2)-unclass(bin1)]
  mat[,decaymat:=mean(decaymat,na.rm = T),by=c("name","diag.idx")]
  mat[,diag.idx:=NULL]
  return(mat)
}

#' Group binned matrices of datasets
#'
#' @param cs CSnorm object, normalized.
#' @param resolution integer. The desired resolution of the matrix.
#' @param group The type of grouping to be performed. Any combination of the given arguments is possible.
#' @param verbose verbose
#' @param ncores integer. The number of cores to parallelize the zeros calculation on.
#'
#' @return CSnorm object
#' @export
#'
#' @examples " "
group_datasets = function(cs, resolution, group=c("condition","replicate","enzyme","experiment"),
                          verbose=T, ncores=1) {
  ### fetch and check inputs
  if (group!="all") group=match.arg(group, several.ok=T)
  if (get_cs_group_idx(cs, resolution=resolution, group=group, raise=F)>0)
    stop("Refusing to overwrite already existing ", group, " grouping.")
  #
  ### compute all matrices
  if (verbose==T) cat("*** compute observed and expected quantities at",resolution,
                      "kb with group=",group,"\n")
  #
  if (verbose==T) cat("   Get zeros per bin\n")
  if (resolution < cs@settings$base.res)
    cat("Warning: resolution is smaller than base.res, results might be inconsistent\n")
  if (resolution == cs@settings$base.res) {
    sbins = cs@settings$sbins
    zeros = cs@zeros
  } else {
    sbins = seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
    zeros = binless:::get_nzeros(cs, sbins, ncores=ncores)
  }
  #
  #predict means, put in triangular form, add biases, and add signal column if absent
  if (verbose==T) cat("   Predict means\n")
  cts.common = binless:::gauss_common_muhat_mean(cs, zeros, sbins)
  cts = binless:::gauss_signal_muhat_mean(cs, cts.common)
  eCmat = cs@design[,.(name,eC=cs@par$eC)]
  cts = merge(cts, eCmat, by="name")
  cts[,log_bias:=lmu.nosig-log_decay-eC]
  #
  if (verbose==T) cat("   Group\n")
  if (group=="all") {
    names=cs@experiments[,unique(name)]
    groups=data.table(name=names,groupname=names)
  } else {
    groups=cs@experiments[,.(name,groupname=do.call(paste,mget(group))),by=group][,.(name,groupname)] #we already know groupname is unique
    groups[,groupname:=ordered(groupname)] #same class as name
  }
  setkey(groups,name)
  cts = groups[cts]
  cts[,name:=groupname]
  cts[,groupname:=NULL]
  setkeyv(cts,c("name","bin1","bin2"))
  #
  if (verbose==T) cat("*** build binned matrices for each experiment\n")
  mat = predict_binned_matrices_irls(copy(cts), cs@par$alpha)
  setkey(mat,name,bin1,bin2)
  #
  if (verbose==T) cat("*** write begin/end positions\n")
  mat = add_bin_bounds_and_distance(mat)
  ### store matrices
  csg=new("CSgroup", mat=mat, interactions=list(), resolution=resolution, group=group,
          cts=cts, par=list(alpha=cs@par$alpha, dmin=cs@settings$dmin, nbins=length(sbins)-1, sbins=sbins),
          names=groups)
  cs@groups=append(cs@groups,list(csg))
  return(cs)
}

#' Bin normalized datasets
#' 
#' @inheritParams group_datasets
#'   
#' @export
#' 
#' @examples " "
bin_all_datasets = function(cs, resolution=cs@settings$base.res, ncores=1, verbose=T) {
  group_datasets(cs, resolution=resolution, group="all", ncores=ncores, verbose=verbose)
}

