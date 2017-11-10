#' @include binless.R
NULL

#' Apply the ICE algorithm to a binned matrix
#' 
#' TODO
#'
#' @param mat a matrix obtained by grouping
#' @param niterations positive integer. Number of iterations to perform
#'
#' @return a CSbinned object containing the ICEd matrix
#' @export
#'
#' @examples
iterative_normalization = function(mat, niterations=100, namecol="name", verbose=T) {
  if (verbose==T) cat("*** iterative normalization with ",niterations," iterations\n")
  raw=mat[,.(name,bin1,bin2,observed)]
  setkey(raw,name,bin1,bin2)
  binned = foreach (n=raw[,unique(get(namecol))], .combine="rbind") %do% {
    binned = raw[get(namecol)==n&bin1<bin2,.(bin1,bin2,N=observed)]
    binned = rbind(binned, binned[,.(bin1=bin2,bin2=bin1,N)])
    binned[,N.weighted:=N]
    #iterate
    for (i in 1:niterations) {
      binned[,b1:=sum(N.weighted),by=bin1]
      binned[,b2:=sum(N.weighted),by=bin2]
      binned[,b1:=b1/mean(b1)]
      binned[,b2:=b2/mean(b2)]
      binned[,N.weighted:=N.weighted/b1/b2]
    }
    binned[,c("b1","b2","N"):=list(NULL,NULL,NULL)]
    setnames(binned,"N.weighted",paste0("ice.",niterations))
    binned=binned[bin1<bin2]
    binned[,c(namecol):=n]
    setkeyv(binned,c(namecol,"bin1","bin2"))
  }
  if ("begin1" %in% names(mat)) 
    binned=merge(binned,mat[,.(name,bin1,begin1,end1,bin2,begin2,end2)],
                 by=c("name","bin1","bin2"),all.x=T)
  setkey(binned,name,bin1,bin2)
  binned
}

#' Perform binning of data
#' 
#' fast IRLS and zero counts approximation
#'
#' @param cts
#' @param dispersion
#' @param verbose boolean.
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
predict_binned_matrices_irls = function(cts, dispersion, verbose=T) {
  #matrices
  if (verbose==T) cat("   Other matrices\n")
  cts[,c("decay","biases","mu.nosig"):=list(exp(log_decay),exp(log_bias),exp(lmu.nosig))]
  mat=cts[,.(ncounts=sum(weight),
             observed=sum(count*weight),
             biasmat=weighted.mean(biases,weight),
             decaymat=weighted.mean(decay,weight),
             background=sum(mu.nosig*weight),
             background.sd=sqrt(sum((mu.nosig+mu.nosig^2/dispersion)*weight)),
             residual=sum(count*weight)/sum(mu*weight),
             normalized=sum(count*weight)/sum(exp(lmu.nosig-log_decay)*weight))
          ,keyby=c("name","bin1","bin2")]
  return(mat)
}

#' Group binned matrices of datasets
#'
#' @param cs CSnorm object, normalized.
#' @param resolution integer. The desired resolution of the matrix.
#' @param group The type of grouping to be performed. Any combination of the given arguments is possible.
#' @param verbose
#' @param ncores integer. The number of cores to parallelize the zeros calculation on.
#'
#' @return CSnorm object
#' @export
#'
#' @examples
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
  mat = predict_binned_matrices_irls(copy(cts), cs@par$alpha, verbose=verbose)
  setkey(mat,name,bin1,bin2)
  #
  if (verbose==T) cat("*** write begin/end positions\n")
  mat = add_bin_begin_and_end(mat)
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
#' @examples
bin_all_datasets = function(cs, resolution=cs@settings$base.res, ncores=1, verbose=T) {
  group_datasets(cs, resolution=resolution, group="all", ncores=ncores, verbose=verbose)
}

