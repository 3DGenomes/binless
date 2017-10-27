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

#' bin matrix and compute data
#' @keywords internal
#' 
gauss_binning_muhat_mean = function(cs, cts.common) {
  cts = cts.common[bin1<=bin2]
  #put in triangular form
  cts2 = cts.common[bin1>bin2]
  setnames(cts2,c("bin1","bin2"),c("bin2","bin1"))
  cts = rbind(cts,cts2)[,.(name,bin1,bin2,dbin,count,lmu.nosig,phi,z,mu,var,log_decay,weight=weight/2)] #each count appears twice
  rm(cts2)
  cts = cts[,.(count=weighted.mean(count,weight/var),lmu.nosig=weighted.mean(lmu.nosig,weight/var),
               z=weighted.mean(z,weight/var),mu=exp(weighted.mean(lmu.nosig+phi,weight/var)),var=1/weighted.mean(1/var,weight),
               log_decay=weighted.mean(log_decay,weight/var),weight=sum(weight)),
            keyby=c("name","bin1","bin2","dbin")]
  stopifnot(cts[,all(bin1<=bin2)])
  return(cts)
}

#' Perform binning of data
#' 
#' fast IRLS and zero counts approximation
#'
#' @param cts
#' @param dispersion
#' @param ncores integer. Number of cores for zero binning
#' @param niter integer. Maximum number of IRLS iterations
#' @param tol numeric. Convergence tolerance for IRLS objective
#' @param verbose boolean.
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
predict_binned_matrices_irls = function(cts, dispersion, ncores=1, niter=100, tol=1e-3, verbose=T) {
  #matrices
  if (verbose==T) cat("   Other matrices\n")
  cts[,c("decay","biases"):=list(exp(log_decay),exp(log_bias))]
  mat=cts[,.(ncounts=sum(weight),
             observed=sum(count*weight),
             biasmat=sum(biases*weight)/sum(weight),
             decaymat=sum(decay*weight)/sum(weight),
             expected=sum(mu*weight),
             expected.sd=sqrt(sum((mu+mu^2/dispersion)*weight)),
             residual=sum(count*weight)/sum(mu*weight))
          ,keyby=c("name","bin1","bin2")]
  #signal matrix
  if (verbose==T) cat("   Signal matrix\n")
  cts[,c("signal","mu.nosig"):=list(1,exp(lmu.nosig))]
  for (i in 1:niter) {
    cts[,c("z","var","signal.old"):=list(count/(signal*mu.nosig)-1,(1/(signal*mu.nosig)+1/dispersion),signal)]
    cts[,signal:=exp(weighted.mean(z+log(signal), weight/var)),by=c("name","bin1","bin2")]
    cts[,signal.sd:=signal[1]*sqrt(1/sum(weight/var)),by=c("name","bin1","bin2")]
    if(cts[,all(abs(signal-signal.old)<tol)]) break
  }
  if (i==niter) cat("Warning: Maximum number of IRLS iterations reached for signal estimation!\n")
  mats = cts[,.(signal=signal[1],signal.sd=signal.sd[1]),keyby=c("name","bin1","bin2")]
  #normalized matrix
  if (verbose==T) cat("   'Normalized' matrix\n")
  mat=mat[mats]
  mat[,normalized:=signal*decaymat]
  mat[,normalized.sd:=signal.sd*decaymat]
  mat[observed==0,c("signal","normalized","signal.sd","normalized.sd"):=list(0,0,NA,NA)]
  return(mat)
}

#' Group binned matrices of datasets
#'
#' @param cs CSnorm object, normalized.
#' @param resolution integer. The desired resolution of the matrix.
#' @param group The type of grouping to be performed. Any combination of the given arguments is possible.
#' @param verbose
#' @param ncores integer. The number of cores to parallelize on.
#'
#' @return CSnorm object
#' @export
#'
#' @examples
group_datasets = function(cs, resolution, group=c("condition","replicate","enzyme","experiment"),
                          verbose=T, ncores=1, niter=100, tol=cs@settings$tol) {
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
  cts = binless:::gauss_binning_muhat_mean(cs, cts.common)
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
  mat = predict_binned_matrices_irls(copy(cts), cs@par$alpha, ncores=ncores, niter=niter, tol=tol, verbose=verbose)
  setkey(mat,name,bin1,bin2)
  #
  if (verbose==T) cat("*** write begin/end positions\n")
  bin1.begin=mat[,bin1]
  bin1.end=mat[,bin1]
  bin2.begin=mat[,bin2]
  bin2.end=mat[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[][,)]")[[3]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[][,)]")[[3]]
  mat[,begin1:=as.integer(as.character(bin1.begin))]
  mat[,end1:=as.integer(as.character(bin1.end))]
  mat[,begin2:=as.integer(as.character(bin2.begin))]
  mat[,end2:=as.integer(as.character(bin2.end))]
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
bin_all_datasets = function(cs, resolution=10000, ncores=1, verbose=T, niter=100, tol=cs@settings$tol) {
  group_datasets(cs, resolution=resolution, group="all", ncores=ncores, verbose=verbose, niter=niter)
}

#' Generate iota and rho genomic biases on evenly spaced points along the genome
#'
#' @param biases data.table.
#' @param beta_iota,beta_rho vectors. spline parameters
#' @param bf_per_kb number of basis functions per kb
#' @param points_per_kb number of evaluation points per kb
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
generate_genomic_biases = function(biases, beta_iota, beta_rho, bf_per_kb=1, points_per_kb=100) {
  begin=biases[,min(pos)]
  end=biases[,max(pos)]
  genome_sz=end-begin
  Krow=round(bf_per_kb*genome_sz/1000)
  S=round(points_per_kb*genome_sz/1000)
  out <- capture.output(
    op<-optimizing(stanmodels$gen_genomic_biases, data=list(Krow=Krow, S=S, begin=begin, end=end,
                                                            beta_iota=beta_iota, beta_rho=beta_rho),
                   as_vector=F, hessian=F, iter=1, verbose=F, init=0))
  dt=as.data.table(op$par)
  setkey(dt, pos)
  return(dt)
}
