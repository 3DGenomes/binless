#' @include csnorm.R
NULL

#' Apply the ICE algorithm to a binned matrix
#'
#' @param csb a CSbinned object
#' @param niterations positive integer. Number of iterations to perform
#'
#' @return a CSbinned object containing the ICEd matrix
#' @export
#'
#' @examples
iterative_normalization = function(raw, niterations=100, namecol="name") {
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
  biases=cs@biases
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


#' Prepare bins and organise into chunks for parallelization
#' @keywords internal
bin_and_chunk = function(cs, resolution, group, ncores) {
  if (group=="all") {
    names=cs@experiments[,unique(name)]
    groups=data.table(name=names,groupname=names)
  } else {
    groups=cs@experiments[,.(name,groupname=do.call(paste,mget(group))),by=group][,.(name,groupname)] #we already know groupname is unique
    groups[,groupname:=ordered(groupname)] #same class as name
  }
  setkey(groups,name)
  #retrieve bin borders
  biases=cs@biases
  counts=cs@counts
  bins=seq(biases[,min(pos)-1],biases[,max(pos)+1+resolution],resolution)
  counts[,c("bin1","bin2"):=list(cut(pos1, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                                 cut(pos2, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12))]
  counts[,c("ibin1","ibin2"):=list(as.integer(bin1)-1,as.integer(bin2)-1)]
  biases[,bin:=cut(pos, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
  biases[,ibin:=as.integer(bin)-1]
  #split across cores
  stepsize=max(2,ceiling(length(bins)/ncores))
  counts[,c("chunk1","chunk2"):=list(ibin1 %/% stepsize, ibin2 %/% stepsize)]
  biases[,chunk:=ibin %/% stepsize]
  chunks=CJ(biases[,(min(chunk):max(chunk))],biases[,(min(chunk):max(chunk))])[V1<=V2]
  return(list(counts=counts,biases=biases,chunks=chunks,groups=groups))
}

#' estimate signal and various matrices for these bins
#' @keywords internal
estimate_signal = function(cs, cts, groups) {
  setkey(cts,name,id1,id2)
  cts=csnorm_predict_all(cs, cts, verbose=F)
  cts=groups[cts]
  cts[,newid:=paste(groupname,ibin1,ibin2)]
  cts=rbind(cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.close,log_mean=log_mean_cclose,log_decay)],
            cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.far,log_mean=log_mean_cfar,log_decay)],
            cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.up,log_mean=log_mean_cup,log_decay)],
            cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.down,log_mean=log_mean_cdown,log_decay)])
  setkey(cts,newid)
  if(cts[,.N]==1) {cbegin=c(1,2)} else {cbegin=c(1,cts[,.(newid,row=.I)][newid!=shift(newid),row],cts[,.N+1])}
  #compute each signal contribution separately
  data=list(G=length(cbegin)-1, N=cts[,.N], cbegin=cbegin, count=cts[,count], log_expected=cts[,log_mean],
            log_decay=cts[,log_decay],alpha=cs@par$alpha)
  output=capture.output(op<-optimizing(csnorm:::stanmodels$predict_binned, data=data, as_vector=F, hessian=T, iter=10000, verbose=F, init=0))
  mat=data.table(name=cts[head(cbegin,-1),groupname], bin1=cts[head(cbegin,-1),bin1], bin2=cts[head(cbegin,-1),bin2],
                 ncounts=op$par$ncounts, observed=op$par$observed, expected=op$par$expected, expected.sd=op$par$expected_sd,
                 decaymat=op$par$decaymat, lpdfr=op$par$lpdfr, lpdfs=op$par$lpdfs, lpdf0=op$par$lpdf0,
                 signal=exp(op$par$log_s), signal.sd=exp(op$par$log_s)*sqrt(as.vector(1/(-head(diag(op$hessian),data$G)))),
                 normalized=exp(op$par$log_r), normalized.sd=exp(op$par$log_r)*sqrt(as.vector(1/(-tail(diag(op$hessian),data$G)))))
  mat[observed==0,c("signal","normalized","signal.sd","normalized.sd"):=list(1,1,0,0)]
  mat
}

#' Perform peak and differential detection through model comparison
#'
#' @param cs csnorm object
#' @param resolution target resolution
#' @param group if grouping is to be performed
#' @param ncores number of cpus to parallelize on
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
csnorm_predict_binned = function(cs, resolution, group, ncores=1) {
  #organise data into bins and chunks
  stuff = csnorm:::bin_and_chunk(cs, resolution, group, ncores)
  chunks=stuff$chunks
  counts=stuff$counts
  biases=stuff$biases
  groups=stuff$groups
  #bin all data and estimate signal
  registerDoParallel(cores=ncores)
  mat = foreach (i=chunks[,V1],j=chunks[,V2], .combine=rbind) %dopar% {
    #get zero-filled portion of counts and predict model on it
    cts=counts[chunk1==i&chunk2==j,.(name,id1,id2,pos1,pos2,contact.close,contact.down,contact.far,contact.up,distance)]
    biases1=biases[chunk==i] #just needed to fill the counts matrix
    biases2=biases[chunk==j]
    if (biases1[,.N]>0 & biases2[,.N]>0) {
      cts=fill_zeros(cts,biases1,biases2,circularize=cs@settings$circularize,dmin=cs@settings$dmin)
      cts=merge(cts,biases1[,.(name,id,pos,bin,ibin)],by.x=c("name","id1","pos1"),by.y=c("name","id","pos"))
      cts=merge(cts,biases2[,.(name,id,pos,bin,ibin)],by.x=c("name","id2","pos2"),by.y=c("name","id","pos"),suffixes=c("1","2"))
      if (cts[,.N]>0) {
        csnorm:::estimate_signal(cs, cts, groups)
      }
    }
  }
  counts[,c("bin1","bin2","ibin1","ibin2"):=list(NULL,NULL,NULL,NULL)]
  biases[,c("bin","ibin"):=list(NULL,NULL)]
  mat
}

#' Predict approximate mean for zero and positive counts at a given binning
#' @keywords internal
csnorm_predict_binned_muhat_irls = function(cs, resolution, zeros) {
  ### predict exact means for positive counts
  init=cs@par
  cpos=copy(cs@counts)
  bsub=cs@biases[,.(name,id,pos)]
  bins=seq(bsub[,min(pos)-1],bsub[,max(pos)+1+resolution],resolution)
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  bsub[,bin:=cut(pos, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
  cpos=merge(bsub[,.(id1=id,bin,log_iota,log_rho)],cpos,by="id1",all.x=F,all.y=T)
  cpos=merge(bsub[,.(id2=id,bin,log_iota,log_rho)],cpos,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  cpos=merge(cbind(cs@design[,.(name)],eC=init$eC), cpos, by="name",all.x=F,all.y=T)
  dbins=cs@settings$dbins
  cpos[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  cpos=merge(cpos,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  cpos[,log_mu.base:=eC + log_decay]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  cpos[,log_mu.base:=NULL]
  #rewrite for each bin
  cpos=rbind(cpos[contact.close>0,.(name, bin1, bin2, cat="close",
                                    count=contact.close, mu=exp(lmu.close), decay=exp(log_decay), weight=1)],
             cpos[contact.far>0,.(name, bin1, bin2, cat="far",
                                  count=contact.far, mu=exp(lmu.far), decay=exp(log_decay), weight=1)],
             cpos[contact.up>0,.(name, bin1, bin2, cat="up",
                                 count=contact.up, mu=exp(lmu.up), decay=exp(log_decay), weight=1)],
             cpos[contact.down>0,.(name, bin1, bin2, cat="down",
                                   count=contact.down, mu=exp(lmu.down), decay=exp(log_decay), weight=1)])
  ### predict approximate means for zero counts
  #approximate decay
  czero = merge(zeros,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  czero = merge(cbind(cs@design[,.(name)],eC=init$eC), czero, by="name",all.x=F,all.y=T)
  czero[,log_mu.base:=eC + log_decay]
  #approximate bias
  bsub = bsub[,.(log_iota=mean(log_iota),log_rho=mean(log_rho)),by=c("name","bin")]
  czero = merge(bsub[,.(name,bin1=bin,log_iota,log_rho)],czero,by=c("name","bin1"),all.x=F,all.y=T)
  czero = merge(bsub[,.(name,bin2=bin,log_iota,log_rho)],czero,by=c("name","bin2"),all.x=F,all.y=T, suffixes=c("2","1"))
  czero[,log_mu:=ifelse(cat=="far", log_mu.base+log_iota1+log_rho2,
                        ifelse(cat=="down", log_mu.base+log_rho1 +log_rho2,
                               ifelse(cat=="close", log_mu.base+log_rho1 +log_iota2,
                                      log_mu.base+log_iota1+log_iota2)))]
  czero[,log_mu.base:=NULL]
  czero = czero[nzero>0,.(name,bin1,bin2,cat,count=0,mu=exp(log_mu),decay=exp(log_decay),weight=nzero)]
  cts=rbind(cpos,czero)
  return(cts)
}

#' Perform peak and differential detection through model comparison
#' 
#' fast IRLS and zero counts approximation
#'
#' @param cs csnorm object
#' @param resolution target resolution
#' @param group if grouping is to be performed
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
csnorm_predict_binned_irls = function(cs, resolution, group, ncores=1, niter=100, tol=1e-3, verbose=T) {
  # get zeros per bin
  if (verbose==T) cat("   Get zeros per bin\n")
  zeros = csnorm:::get_nzeros_binning(cs, resolution, ncores=ncores)
  # predict means
  if (verbose==T) cat("   Predict means\n")
  cts = csnorm:::csnorm_predict_binned_muhat_irls(cs, resolution, zeros)
  # group
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
  #matrices
  if (verbose==T) cat("   Other matrices\n")
  mat=cts[,.(ncounts=sum(weight),
             observed=sum(count*weight),
             expected=sum(mu*weight),
             expected.sd=sqrt(sum((mu+mu^2/init$alpha)*weight)),
             decaymat=sum(decay*weight)/sum(weight),
             lpdf0=sum(dnbinom(count,mu=mu, size=init$alpha, log=T)*weight))
          ,keyby=c("name","bin1","bin2")]
  #signal matrix
  if (verbose==T) cat("   Signal matrix\n")
  cts[,signal:=1]
  for (i in 1:niter) {
    cts[,c("z","var","signal.old"):=list(count/(signal*mu)-1,(1/(signal*mu)+1/init$alpha),signal)]
    cts[,signal:=exp(weighted.mean(z+log(signal), weight/var)),by=c("name","bin1","bin2")]
    cts[,signal.sd:=signal[1]*sqrt(1/sum(weight/var)),by=c("name","bin1","bin2")]
    if(cts[,all(abs(signal-signal.old)<tol)]) break
  }
  if (i==niter) cat("Warning: Maximum number of IRLS iterations reached for signal estimation!\n")
  mats = cts[,.(signal=signal[1],signal.sd=signal.sd[1],
                lpdfs=sum(dnbinom(count,mu=mu*signal, size=init$alpha, log=T)*weight))
             ,keyby=c("name","bin1","bin2")]
  #normalized matrix
  if (verbose==T) cat("   'Normalized' matrix\n")
  cts[,normalized:=1]
  for (i in 1:niter) {
    cts[,c("z","var","normalized.old"):=list(count/(normalized*mu/decay)-1,
                                             (1/(normalized*mu/decay)+1/init$alpha),normalized)]
    cts[,normalized:=exp(weighted.mean(z+log(normalized), weight/var)),by=c("name","bin1","bin2")]
    cts[,normalized.sd:=normalized[1]*sqrt(1/sum(weight/var)),by=c("name","bin1","bin2")]
    if(cts[,all(abs(normalized-normalized.old)<tol)]) break
  }
  if (i==niter) cat("Warning: Maximum number of IRLS iterations reached for normalized estimation!\n")
  matr = cts[,.(normalized=normalized[1],normalized.sd=normalized.sd[1],
                lpdfr=sum(dnbinom(count,mu=mu*normalized/decay, size=init$alpha, log=T)*weight))
             ,keyby=c("name","bin1","bin2")]
  mat=mat[mats[matr]]
  mat[observed==0,c("signal","normalized","signal.sd","normalized.sd"):=list(0,0,NA,NA)]
  return(mat)
}

#' Common call for binning
#' @keywords internal
#'
#' @examples
compute_grouped_matrices = function(cs, resolution, group, ncores, ice, verbose, niter=100, tol=1e-3) {
  if (verbose==T) cat("*** build binned matrices for each experiment\n")
  #mat=csnorm_predict_binned(cs, resolution, group=group, ncores=ncores)
  mat=csnorm_predict_binned_irls(cs, resolution, group=group, ncores=ncores, niter=niter, tol=tol, verbose=verbose)
  setkey(mat,name,bin1,bin2)
  if (ice>0) {
    if (verbose==T) cat("*** iterative normalization with ",ice," iterations\n")
    raw=mat[,.(name,bin1,bin2,observed)]
    setkey(raw,name,bin1,bin2)
    iced=iterative_normalization(raw, niterations=ice, namecol="name")
    setkey(iced,name,bin1,bin2)
    mat=merge(mat,iced,all.x=T,all.y=F)
  }
  #write begins/ends
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
  return(mat)
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

#' Bin normalized datasets
#' 
#' @param cs CSnorm object, optimized.
#' @param resolution integer. The desired resolution of the matrix.
#' @param ncores integer. The number of cores to parallelize on.
#' @param ice integer. If positive, perform the optional Iterative Correction 
#'   algorithm, useful for comparison purposes. The value determines the number
#'   of iterations.
#' @param verbose
#'   
#' @return A CSnorm object containing an additional CSbinned object in cs@binned
#' @export
#' 
#' @examples
bin_all_datasets = function(cs, resolution=10000, ncores=1, ice=-1, verbose=T) {
  if (get_cs_binned_idx(cs,resolution,raise=F)>0) {
    stop("Refusing to overwrite already existing matrices at ", resolution/1000,
         "kb. Use them or remove them from the cs@binned list")
  }
  mat = csnorm:::compute_grouped_matrices(cs, resolution=resolution, group="all", ncores=ncores, ice=ice, verbose=verbose)
  #create CSmatrix and CSbinned object
  csm=new("CSmatrix", mat=mat, group="all", ice=(ice>0), ice.iterations=ice,
          names=as.character(mat[,unique(name)]))
  csb=new("CSbinned", resolution=resolution, grouped=list(csm),
          individual=copy(mat))
  cs@binned=append(cs@binned,csb)
  cs
}

#' Group binned matrices of datasets
#'
#' @param cs CSnorm object
#' @param resolution see \code{\link{bin_all_datasets}}, used to identify the input matrices.
#' @param group The type of grouping to be performed. Any combination of the given arguments is possible.
#' @inheritParams bin_all_datasets
#'
#' @return CSnorm object
#' @export
#'
#' @examples
group_datasets = function(cs, resolution, group=c("condition","replicate","enzyme","experiment"),
                          ice=-1, verbose=T, ncores=1) {
  #fetch and check inputs
  experiments=cs@experiments
  csbi=get_cs_binned_idx(cs, resolution=resolution, raise=T)
  csb=cs@binned[[csbi]]
  group=match.arg(group, several.ok=T)
  if (get_cs_matrix_idx(csb, group, raise=F)>0) {
    stop("Refusing to overwrite already existing ", group,
         " group matrices. Use them or remove them from the cs@binned[[",csbi,
         "]]@grouped list and @metadata table")
  }
  #
  mat = compute_grouped_matrices(cs, resolution=resolution, group=group, ncores=ncores, ice=ice, verbose=verbose)
  #store matrices
  csm=new("CSmatrix", mat=mat, group=group, ice=(ice>0), ice.iterations=ice,
          names=as.character(mat[,unique(name)]))
  csb@grouped=append(csb@grouped,list(csm))
  cs@binned[[csbi]]=csb
  return(cs)
}
