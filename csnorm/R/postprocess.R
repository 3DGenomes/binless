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

#' Prepare bins and organise into chunks for parallelization
#' @keywords internal
bin_and_chunk = function(cs, resolution, group, ncores) {
  if (group=="all") {
    names=cs@experiments[,unique(name)]
    groups=data.table(name=names,groupname=names)
  } else {
    groups=cs@experiments[,.(name,groupname=do.call(paste,mget(group))),by=group][,.(name,groupname)] #we already know groupname is unique
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
  stepsize=max(2,ceiling(length(bins)/(5*ncores)))
  counts[,c("chunk1","chunk2"):=list(ibin1 %/% stepsize, ibin2 %/% stepsize)]
  biases[,chunk:=ibin %/% stepsize]
  chunks=CJ(biases[,min(chunk):max(chunk)],biases[,min(chunk):max(chunk)])[V1<=V2]
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
    cts=counts[chunk1==i&chunk2==j]
    biases1=biases[chunk==i] #just needed to fill the counts matrix
    biases2=biases[chunk==j]
    if (biases1[,.N]>0 & biases2[,.N]>0) {
      cts=fill_zeros(cts,biases1,biases2,circularize=cs@settings$circularize,dmin=cs@settings$dmin)
      if (cts[,.N]>0) {
        csnorm:::estimate_signal(cs, cts, groups)
      }
    }
  }
  counts[,c("bin1","bin2","ibin1","ibin2"):=list(NULL,NULL,NULL,NULL)]
  biases[,c("bin","ibin"):=list(NULL,NULL)]
  mat
}

#' Common call for binning
#' @keywords internal
#'
#' @examples
bin_grouped_matrix = function(cs, resolution, group, ncores, ice, verbose) {
  if (verbose==T) cat("*** build binned matrices for each experiment\n")
  mat=csnorm_predict_binned(cs, resolution, group="all", ncores=ncores)
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

#' compute log p(D|s)
#' @keywords internal
get_signal_likelihood = function(mat, chunks, counts, biases, groups, ncores) {
  #bin all data and estimate signal
  registerDoParallel(cores=ncores)
  liks = foreach (i=chunks[,V1],j=chunks[,V2], .combine=rbind) %dopar% {
    #get zero-filled portion of counts and predict model on it
    cts=counts[chunk1==i&chunk2==j]
    biases1=biases[chunk==i] #just needed to fill the counts matrix
    biases2=biases[chunk==j]
    if (biases1[,.N]>0 & biases2[,.N]>0) {
      cts=fill_zeros(cts,biases1,biases2,circularize=cs@settings$circularize,dmin=cs@settings$dmin)
      if (cts[,.N]>0) {
        setkey(cts,name,id1,id2)
        cts=csnorm_predict_all(cs, cts, verbose=F)
        cts=groups[cts]
        cts[,newid:=paste(groupname,ibin1,ibin2)]
        cts=rbind(cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.close,log_mean=log_mean_cclose)],
                  cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.far,log_mean=log_mean_cfar)],
                  cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.up,log_mean=log_mean_cup)],
                  cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.down,log_mean=log_mean_cdown)])
        setkey(cts,newid)
        signals=merge(mat,cts[,.SD[1],by=newid], by.x=c("name","bin1","bin2"),by.y=c("groupname","bin1","bin2"), sort=F)
        setkey(signals,newid)
        if(cts[,.N]==1) {cbegin=c(1,2)} else {cbegin=c(1,cts[,.(newid,row=.I)][newid!=shift(newid),row],cts[,.N+1])}
        #compute each signal contribution separately
        data=list(G=length(cbegin)-1, N=cts[,.N], cbegin=cbegin, count=cts[,count], log_expected=cts[,log_mean],
                  alpha=cs@par$alpha, log_s=signals[,log_s], S=10)
        output=capture.output(op<-optimizing(csnorm:::stanmodels$likelihood_binned, data=data, as_vector=F, hessian=F, iter=0, verbose=F, init=0))
        out = foreach (g=1:data$G,.combine=rbind) %do%
          data.table(bin=g,logp=op$par$lpdfs[g,],log_sref=data$log_s[g], lpdf0=op$par$lpdf0[g], lpdfref=op$par$lpdfref[g], idx=1:data$S)
        out[,log_s:=log_sref*(idx-1)/(data$S-1)]
        out[,expected:=(lpdf0-lpdfref)*(log_s/log_sref-1)^2+lpdfref]
        ggplot(out)+geom_point(aes(log_s,logp))+facet_wrap(~bin,scales="free")+
          geom_line(aes(log_s,expected),colour="red")
        
        lik=data.table(name=cts[head(cbegin,-1),groupname], bin1=cts[head(cbegin,-1),bin1], bin2=cts[head(cbegin,-1),bin2],
                       ncounts=op$par$ncounts, observed=op$par$observed, expected=op$par$expected, expected.sd=op$par$expected_sd,
                       decaymat=op$par$decaymat, lpdfr=op$par$lpdfr, lpdfs=op$par$lpdfs, lpdf0=op$par$lpdf0,
                       signal=exp(op$par$log_s), signal.sd=exp(op$par$log_s)*sqrt(as.vector(1/(-head(diag(op$hessian),data$G)))),
                       normalized=exp(op$par$log_r), normalized.sd=exp(op$par$log_r)*sqrt(as.vector(1/(-tail(diag(op$hessian),data$G)))))
        lik[observed==0,c("signal","normalized","signal.sd","normalized.sd"):=list(0,0,0,0)]
        lik
      }
    }
  }
}

#' Perform peak calling through lasso and model comparison
#' @keywords internal
interactions_binned = function(cs, resolution=resolution, group=group, threshold=threshold,
                               ncores=ncores, alpha.gamma=2, beta.gamma=1) {
  #initial gamma and log_s estimates
  likmat = get_matrices(cs, resolution=resolution, group=group)[bin2>=bin1,.(name,bin1,bin2,log_sref=log(signal),
                                                                             lpdf0,lpdfref=lpdfs, signal.sd)]
  #gamma = likmat[,.N/sum(abs(log_sref))]
  likmat[,log_s:=log_sref]
  nbins=likmat[,.N]
  #alternate recalculation of log-likelihood and estimation of gamma
  registerDoParallel(cores=ncores)
  bfac = foreach (gamma=likmat[,unique(abs(log_sref))], .combine=rbind) %dopar% {
    #soft-thresholding of log_s
    mat=copy(likmat)
    mat[,log_s:=sign(log_sref)*pmax(abs(log_sref)-gamma,0)]
    #interpolate log likelihood using laplace approximation
    mat[,lpdfs:=ifelse(log_sref==0,lpdf0,(lpdf0-lpdfref)*(log_s/log_sref-1)^2+lpdfref)]
    #compute K(gamma)
    data.table(gamma=gamma,nnz=mat[abs(log_sref)>gamma,.N],params=paste(alpha.gamma,beta.gamma),
               logK=mat[,sum(lpdfs-lpdf0-gamma*(abs(log_s)+beta.gamma)+alpha.gamma*log(beta.gamma*gamma)
                             - log(2) - lgamma(alpha.gamma))])
  }
  #select best model and produce matrix with interactions
  bfac=bfac[logK==max(logK)]
  retlist=list(gamma=bfac[,gamma], prob.signal=bfac[,1/(exp(-logK)+1)])
  retlist$is.significant=retlist$prob.signal>threshold
  if(retlist$is.significant==TRUE) {
    mat=likmat[,.(name,bin1,bin2,signal=exp(sign(log_sref)*pmax(abs(log_sref)-retlist$gamma,0)),signal.sd)]
  } else {
    mat=likmat[,.(name,bin1,bin2,signal=1,signal.sd=0)]
  }
  mat[,direction:=ifelse(signal>1,"enriched",ifelse(signal==1,NA,"depleted"))]
  retlist$mat=mat
  retlist
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
  mat = bin_grouped_matrix(cs, resolution=resolution, group="all", ncores=ncores, ice=ice, verbose=verbose)
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
group_datasets = function(cs, resolution, group=c("condition","replicate","enzyme","experiment"), ice=-1, verbose=T) {
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
  mat = bin_grouped_matrix(cs, resolution=resolution, group=group, ncores=ncores, ice=ice, verbose=verbose)
  #store matrices
  csm=new("CSmatrix", mat=mat, group=group, ice=(ice>0), ice.iterations=ice,
          names=as.character(mat[,unique(name)]))
  csb@grouped=append(csb@grouped,list(csm))
  cs@binned[[csbi]]=csb
  return(cs)
}

#' Detect significant interactions wrt expected
#' 
#' @param cs CSnorm object
#' @param resolution,group see
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to
#'   identify the input matrices.
#' @param threshold significance threshold, between 0 and 1
#' @param ncores number of cores used for parallelization
#' @inheritParams call_interactions
#'   
#' @return the binned matrix with additional information relating to these 
#'   significant interactions
#' @export
#' 
#' @examples
detect_interactions = function(cs, resolution, group, threshold=0.95, ncores=1){
  #get CSmat object
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csm, type="interactions", threshold=threshold, ref="expected", raise=F)>0) {
    stop("Refusing to overwrite this already detected interaction")
  }
  stuff = interactions_binned(cs, resolution=resolution, group=group, threshold=threshold, ncores=ncores)
  mat=stuff$mat
  #store back
  csi=new("CSinter", mat=stuff$mat, type="interactions", threshold=threshold, ref="expected",
          gamma=stuff$gamma, prob.signal=stuff$prob.signal, nnz=mat[signal!=1,.N])
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}


#' Detect significant differences with a reference
#' 
#' @param binned as returned by \code{\link{csnorm_predict_binned}}
#' @inheritParams call_interactions
#'   
#' @return the binned matrix with additional information relating to these
#'   significant interactions
#' @export
#' 
#' @examples
detect_differences = function(cs, resolution, group, ref, threshold=0.95, ncores=1){
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csm, type="differences", threshold=threshold, ref="expected", raise=F)>0) {
    stop("Refusing to overwrite this already detected interaction")
  }
  mat=detection_binned(cs, resolution=resolution, group=group, ref=ref, threshold=threshold, ncores=ncores)
  #add interaction to cs object
  csi=new("CSinter", mat=mat, type="differences", threshold=threshold, ref=ref)
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}
  
