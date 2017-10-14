#' @include csnorm.R
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

#' count number of zeros in each rectangular bin
#' @keywords internal
#' 
get_nzeros_binning = function(cs, resolution, ncores=1) {
  stopifnot(cs@counts[id1>=id2,.N]==0)
  #count left and right
  cts=melt(cs@counts[,.(name,pos1,pos2,distance,contact.close,contact.down,contact.far,contact.up)],
           id.vars=c("name","pos1","pos2","distance"))[value>0]
  #retrieve bin borders
  biases=cs@biases[,.(name,id,pos)]
  bins=seq(biases[,min(pos)-1],biases[,max(pos)+1+resolution],resolution)
  biases[,bin:=cut(pos, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
  cts[,c("bin1","bin2","dbin"):=
        list(cut(pos1, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
             cut(pos2, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
             cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  #count per bin
  cts=cts[,.(nnz=.N),keyby=c("name","bin1","bin2","dbin","variable")]
  #Count the number of crossings per distance bin
  #looping over IDs avoids building NxN matrix
  registerDoParallel(cores=ncores)
  chunksize=cs@biases[,ceiling(.N/(10*ncores))]
  nchunks=cs@biases[,ceiling(.N/chunksize)]
  crossings = foreach(chunk=1:nchunks, .combine=rbind) %dopar% {
    bs=biases[((chunk-1)*chunksize+1):min(.N,chunk*chunksize)]
    foreach(n=bs[,name], p=bs[,pos], b=bs[,bin], .combine=rbind) %do% {
      crossings = biases[name==n&pos>p,.(name,bin2=bin,distance=abs(pos-p))]
      if (cs@settings$circularize>0)  crossings[,distance:=pmin(distance,cs@settings$circularize+1-distance)]
      crossings[,dbin:=cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
      crossings[distance>=cs@settings$dmin,.(bin1=b,ncross=.N),by=c("name","bin2","dbin")]
    }
  }
  stopImplicitCluster()
  crossings=crossings[,.(ncross=sum(ncross)),keyby=c("name","bin1","bin2","dbin")]
  zeros = rbind(
    merge(crossings,cts[variable=="contact.close"],by=c("name","bin1","bin2","dbin"),all=T)[
    ,.(name,bin1,bin2,dbin,cat="close",ncross,nnz)],
    merge(crossings,cts[variable=="contact.far"],by=c("name","bin1","bin2","dbin"),all=T)[
      ,.(name,bin1,bin2,dbin,cat="far",ncross,nnz)],
    merge(crossings,cts[variable=="contact.down"],by=c("name","bin1","bin2","dbin"),all=T)[
      ,.(name,bin1,bin2,dbin,cat="down",ncross,nnz)],
    merge(crossings,cts[variable=="contact.up"],by=c("name","bin1","bin2","dbin"),all=T)[
      ,.(name,bin1,bin2,dbin,cat="up",ncross,nnz)])
  zeros[is.na(nnz),nnz:=0]
  zeros[,nzero:=ncross-nnz]
  stopifnot(zeros[is.na(ncross),.N==0])
  stopifnot(zeros[nzero<0,.N==0])
  return(zeros)
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
csnorm_predict_binned_matrices_irls = function(cts, dispersion, ncores=1, niter=100, tol=1e-3, verbose=T) {
  #matrices
  if (verbose==T) cat("   Other matrices\n")
  cts[,c("decay","biases"):=list(exp(log_decay),exp(log_bias))]
  mat=cts[,.(ncounts=sum(weight),
             observed=sum(count*weight),
             biasmat=sum(biases*weight)/sum(weight),
             decaymat=sum(decay*weight)/sum(weight),
             expected=sum(mu*weight),
             expected.sd=sqrt(sum((mu+mu^2/dispersion)*weight)),
             residual=sum(count*weight)/sum(mu*weight),
             lpdf0=sum(dnbinom(count,mu=mu, size=dispersion, log=T)*weight))
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
  mats = cts[,.(signal=signal[1],signal.sd=signal.sd[1],
                lpdfs=sum(dnbinom(count,mu=mu.nosig*signal, size=dispersion, log=T)*weight))
             ,keyby=c("name","bin1","bin2")]
  #normalized matrix
  if (verbose==T) cat("   'Normalized' matrix\n")
  cts[,normalized:=decay]
  for (i in 1:niter) {
    cts[,c("z","var","normalized.old"):=list(count/(normalized*mu.nosig/decay)-1,
                                             (1/(normalized*mu.nosig/decay)+1/dispersion),normalized)]
    cts[,normalized:=exp(weighted.mean(z+log(normalized), weight/var)),by=c("name","bin1","bin2")]
    cts[,normalized.sd:=normalized[1]*sqrt(1/sum(weight/var)),by=c("name","bin1","bin2")]
    if(cts[,all(abs(normalized-normalized.old)<tol)]) break
  }
  if (i==niter) cat("Warning: Maximum number of IRLS iterations reached for normalized estimation!\n")
  matr = cts[,.(normalized=normalized[1],normalized.sd=normalized.sd[1],
                lpdfr=sum(dnbinom(count,mu=mu.nosig*normalized/decay, size=dispersion, log=T)*weight))
             ,keyby=c("name","bin1","bin2")]
  mat=mat[mats[matr]]
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
                          verbose=T, ncores=1, niter=100, tol=1e-3) {
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
    zeros = csnorm:::get_nzeros(cs, sbins, ncores=ncores)
  }
  #
  #predict means, put in triangular form, add biases, and add signal column if absent
  if (verbose==T) cat("   Predict means\n")
  cts.common = csnorm:::csnorm_gauss_common_muhat_mean(cs, zeros, sbins)
  cts = csnorm:::csnorm_gauss_signal_muhat_mean(cs, cts.common, zeros, sbins)
  eCmat = cs@design[,.(name,eC=cs@par$eC)]
  cts = merge(cts, eCmat, by="name")
  cts[,log_bias:=lmu.nosig-log_decay-eC]
  if (!("phi" %in% names(cts))) cts[,phi:=0]
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
  mat = csnorm_predict_binned_matrices_irls(copy(cts), cs@par$alpha, ncores=ncores, niter=niter, tol=tol, verbose=verbose)
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
bin_all_datasets = function(cs, resolution=10000, ncores=1, verbose=T, niter=100) {
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
