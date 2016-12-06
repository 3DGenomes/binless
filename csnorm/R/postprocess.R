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

#' Common call for binning
#' @keywords internal
#'
#' @examples
compute_grouped_matrices = function(cs, resolution, group, ncores, ice, verbose) {
  if (verbose==T) cat("*** build binned matrices for each experiment\n")
  mat=csnorm_predict_binned(cs, resolution, group=group, ncores=ncores)
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

#' common calculation between interactions and differences
#' @keywords internal
estimate_significant_common_binned = function(cs, cts, groups, prior.sd) {
  setkey(cts,name,id1,id2)
  cts=csnorm_predict_all(cs, cts, verbose=F)
  cts=groups[cts]
  cts[,newid:=paste(groupname,ibin1,ibin2)]
  cts=rbind(cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.close,log_mean=log_mean_cclose)],
            cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.far,log_mean=log_mean_cfar)],
            cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.up,log_mean=log_mean_cup)],
            cts[,.(newid,groupname,ibin1,ibin2,bin1,bin2,count=contact.down,log_mean=log_mean_cdown)])
  setkey(cts,newid)
  if(cts[,.N]==1) {cbegin=c(1,2)} else {cbegin=c(1,cts[,.(newid,row=.I)][newid!=shift(newid),row],cts[,.N+1])}
  #compute each signal contribution separately
  data=list(G=length(cbegin)-1, N=cts[,.N], cbegin=cbegin, observed=cts[,count], log_expected=cts[,log_mean],
            alpha=cs@par$alpha, sigma=prior.sd)
  output=capture.output(op<-optimizing(csnorm:::stanmodels$detection, data=data, as_vector=F, hessian=T, iter=1000, verbose=F, init=0))
  return(list(cts=cts,op=op,cbegin=cbegin))
}

#' estimate significant interactions for a set of counts
#' @keywords internal
estimate_significant_interactions_binned = function(cs, cts, groups, threshold, prior.sd) {
  stuff=estimate_significant_common_binned(cs,cts,groups, prior.sd)
  op=stuff$op
  cts=stuff$cts
  cbegin=stuff$cbegin
  mat=data.table(name=cts[head(cbegin,-1),groupname], bin1=cts[head(cbegin,-1),bin1], bin2=cts[head(cbegin,-1),bin2],
                 K=exp(op$par$lpdfs+1/2*log(2*pi*as.vector(1/(-diag(op$hessian))))-op$par$lpdf0),
                 direction=ifelse(0<op$par$log_s,"enriched","depleted"),
                 signal.signif=exp(op$par$log_s))
  mat[,signal.signif.sd:=sqrt(as.vector(1/(-diag(op$hessian))))*signal.signif]
  mat[op$par$binned==0,c("signal.signif","signal.signif.sd","K","direction"):=list(1,0,0,NA)]
  mat[,prob.gt.expected:=K/(1+K)]
  mat[,K:=NULL]
  mat[,is.significant:=prob.gt.expected > threshold]
  mat
}

#' estimate significant differences wrt ref for a set of counts
#' @keywords internal
estimate_significant_differences_binned = function(cs, cts, groups, ref, threshold, prior.sd) {
  stuff=estimate_significant_common_binned(cs,cts,groups, prior.sd)
  op=stuff$op
  cts=stuff$cts
  cbegin=stuff$cbegin
  #compute grouped signal contributions
  mat=data.table(name=cts[head(cbegin,-1),groupname], bin1=cts[head(cbegin,-1),bin1], bin2=cts[head(cbegin,-1),bin2],
                 lpdms=op$par$lpdfs+1/2*log(2*pi*as.vector(1/(-diag(op$hessian)))), #log(p(D|Msignal))
                 log_s.mean=op$par$log_s, signal=exp(op$par$log_s))
  mat[,signal.sd:=sqrt(as.vector(1/(-diag(op$hessian))))*signal]
  mat[op$par$binned==0,c("signal","signal.sd","lpdms","log_s.mean"):=list(1,0,op$par$lpdf0[op$par$binned==0],0)]
  refcounts=cts[groupname==ref]
  diffcounts=foreach(n=cts[groupname!=ref,unique(groupname)],.combine=rbind) %do% {
    tmp=rbind(cts[groupname==n],refcounts)
    tmp[,groupname:=n]
    tmp[,newid:=paste(n,ibin1,ibin2)]
  }
  setkey(diffcounts,newid)
  if(diffcounts[,.N]==1){cbegin=c(1,2)} else {cbegin=c(1,diffcounts[,.(newid,row=.I)][newid!=shift(newid),row],diffcounts[,.N+1])}
  data=list(G=length(cbegin)-1, N=diffcounts[,.N], cbegin=cbegin,
            observed=diffcounts[,count], log_expected=diffcounts[,log_mean], alpha=cs@par$alpha, sigma=prior.sd)
  output=capture.output(op2<-optimizing(csnorm:::stanmodels$detection, data=data, as_vector=F, hessian=T, iter=1000, verbose=F, init=0))
  refmat=mat[name==ref,.(bin1,bin2,lpdmref=lpdms,logsref=log_s.mean,refsignal=signal,refsignal.sd=signal.sd)]
  mat=merge(mat[name!=ref],refmat,by=c("bin1","bin2"))
  mat[,c("lpdm2","direction"):=list(lpdms+lpdmref,ifelse(logsref<log_s.mean,"enriched","depleted"))]
  mat[,c("lpdms","lpdmref","logsref","log_s.mean"):=list(NULL,NULL,NULL,NULL)]
  mat2=data.table(name=diffcounts[head(cbegin,-1),groupname],
                  bin1=diffcounts[head(cbegin,-1),bin1], bin2=diffcounts[head(cbegin,-1),bin2],
                  lpdms=op2$par$lpdfs+1/2*log(2*pi*as.vector(1/(-diag(op2$hessian)))))
  mat2[op2$par$binned==0,lpdms:=op2$par$lpdf0]
  mat=merge(mat,mat2,by=c("name","bin1","bin2"))
  mat[,K:=exp(lpdm2-lpdms)]
  colname=paste0("prob.gt.",ref)
  mat[,c(colname):=K/(1+K)]
  mat[,is.significant:=get(colname) > threshold]
  mat[,c("difference","difference.sd"):=list(signal/refsignal,(signal.sd*refsignal-refsignal.sd*signal)/(signal*refsignal))]
  mat[,c("lpdm2","lpdms","K","signal","signal.sd","refsignal","refsignal.sd"):=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL)]
  mat
}

#' Perform peak and differential detection through model comparison
#'
#' @param cs 
#' @param ref 
#' @param resolution 
#' @param group 
#' @param threshold on the probability K/(1+K) where K is the Bayes factor
#' @param ncores 
#' @param prior.sd 
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
csnorm_detect_binned = function(cs, resolution, group, ref="expected", threshold=0.95, prior.sd=5, ncores=1) {
  stuff = bin_and_chunk(cs, resolution, group, ncores)
  chunks=stuff$chunks
  counts=stuff$counts
  biases=stuff$biases
  groups=stuff$groups
  if (ref != "expected") {
    if (!(ref %in% groups[,groupname]))
      stop("Reference group name not found! Valid names are: ",deparse(as.character(groups[,groupname])))
    if (groups[groupname!=ref,.N]==0)
      stop("There is no other group than ",ref, ", cannot compute differences!")
  }
  #bin all data and 
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
        if (ref=="expected") {
          csnorm:::estimate_significant_interactions_binned(cs,cts,groups, threshold, prior.sd)
        } else {
          csnorm:::estimate_significant_differences_binned(cs,cts,groups,ref, threshold, prior.sd)
        }
      }
    }
  }
  counts[,c("bin1","bin2","ibin1","ibin2"):=list(NULL,NULL,NULL,NULL)]
  biases[,c("bin","ibin"):=list(NULL,NULL)]
  mat
}

#' group counts to compute signal by cross-validated lasso
#' @keywords internal
estimate_binless_signal = function(cs, cts, groups, mat) {
  setkey(cts,name,id1,id2)
  cts=csnorm_predict_all(cs, cts, verbose=F)
  cts=groups[cts]
  cts=rbind(cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.close,log_mean=log_mean_cclose)],
            cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.far,log_mean=log_mean_cfar)],
            cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.up,log_mean=log_mean_cup)],
            cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.down,log_mean=log_mean_cdown)])
  cts=mat[,.(groupname,ibin1,ibin2,phi)][cts,,on=c("groupname","ibin1","ibin2")]
  cts[,var:=1/exp(log_mean)+1/cs@par$alpha]
  cts[,phihat:=(count/exp(log_mean)-1+phi)]
  #compute on all counts
  ret=cts[,.(phihat=weighted.mean(phihat,1/var),var=1/sum(1/var),ncounts=.N),
          by=c("groupname","ibin1","ibin2","bin1","bin2","phi")]
  return(ret)
}

#' connectivity on a triangle
#' @keywords internal
flsa_compute_connectivity = function(nbins, start=0) {
  stopifnot(start>=0,nbins>=2)
  if (nbins==2) {
    ret=sapply(list(start+1,c(start,start+2),start+1),as.integer)
    names(ret) = c(start, start+1, start+2)
  } else {
    upper.row = c(list(start+1),
                  sapply((start+1):(start+nbins-2),function(x){c(x-1,x+1,x+nbins-1)},simplify=F),
                  list(c(start+nbins-2,start+2*(nbins-1))))
    names(upper.row) = start:(start+nbins-1)
    lower.tri = flsa_compute_connectivity(nbins-1,start=start+nbins)
    for (i in 1:(nbins-1))
      lower.tri[[i]] = c(lower.tri[[i]], start+i)
    ret = sapply(c(upper.row,lower.tri),as.integer)
  }
  class(ret) = "connListObj"
  return(ret)
}

#' return phi coefficients for a given value of lambda
#' @keywords internal
flsa_get_phi = function(res,lambda,mat) {
  mat[,value:=flsaGetSolution(res, lambda1=0, lambda2=lambda)[1,]]
  mat[,phi:=value*sqrt(var)]
  return(sol)
  #plot for a specific value of lambda
  ggplot()+geom_raster(aes(ibin1,ibin2,fill=value),data=mat)
  ggplot()+geom_raster(aes(ibin1,ibin2,fill=phi),data=mat)+
    geom_raster(aes(ibin2,ibin1,fill=phihat),data=mat[groupname==groupname[1]])
}

#' Perform binless interaction detection using fused lasso
#'
#' @param cs 
#' @param ref 
#' @param resolution 
#' @param group 
#' @param threshold on the probability K/(1+K) where K is the Bayes factor
#' @param ncores 
#' @param niter number of IRLS iterations
#' @param cv.fold perform x-fold cross-validation
#' @param cv.gridsize compute possible values of log lambda on a grid
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
csnorm_detect_binless = function(cs, resolution, group, ref="expected", threshold=0.95, niter=1,
                                 ncores=1, cv.fold=10, cv.gridsize=100) {
  stuff = csnorm:::bin_and_chunk(cs, resolution, group, ncores)
  chunks=stuff$chunks
  counts=stuff$counts
  biases=stuff$biases
  groups=stuff$groups
  if (ref != "expected") {
    if (!(ref %in% groups[,groupname]))
      stop("Reference group name not found! Valid names are: ",deparse(as.character(groups[,groupname])))
    if (groups[groupname!=ref,.N]==0)
      stop("There is no other group than ",ref, ", cannot compute differences!")
  }
  #fused lasso loop
  registerDoParallel(cores=ncores)
  mat = biases[,CJ(groupname=groups[,unique(groupname)],
                   ibin1=min(ibin):max(ibin),ibin2=min(ibin):max(ibin))][ibin1<=ibin2]
  mat[,phi:=0]
  for (step in 1:niter) {
    #compute signal matrix and new weights in parallel
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
          stopifnot(ref=="expected")
          csnorm:::estimate_binless_signal(cs, cts, groups, mat)
        }
      }
    }
    #split into cross-validation groups and sort
    mat[,cv.group:=as.character(((ibin2+ibin1*max(ibin1))%%cv.fold)+1)]
    setkey(mat,groupname,ibin1,ibin2,bin1,bin2)
    
    #ggplot(mat[groupname=="GM MboI 1"])+geom_raster(aes(ibin1,ibin2,fill=phihat))
    #
    #perform one fused lasso per dataset and cv group
    cvgroups=as.character(1:cv.fold)
    groupnames=groups[,unique(groupname)]
    cmat = csnorm:::flsa_compute_connectivity(mat[,max(ibin2)+1])
    res.ref = foreach (g=groupnames, .final=function(x){setNames(x,groupnames)}) %dopar% {
      submat=mat[groupname==g]
      stopifnot(length(cmat)==submat[,.N]) #ibin indices have something wrong
      flsa(submat[,phihat/sqrt(var)], connListObj=cmat, verbose=F)
    }
    res.cv = foreach (cv=cvgroups, .final=function(x){setNames(x,cvgroups)}) %:%
      foreach (g=groupnames, .final=function(x){setNames(x,groupnames)}) %dopar% {
        submat=copy(mat[groupname==g])
        submat[cv.group==cv,phihat:=0]
        stopifnot(length(cmat)==submat[,.N]) #ibin indices have something wrong
        flsa(submat[,phihat/sqrt(var)], connListObj=cmat, verbose=F)
    }
    #get grid values for lambda
    lambdas = unique(sort(c(sapply(res.ref,function(x){x$BeginLambda}),recursive=T)))
    lambdas=lambdas[lambdas>7]
    lmin = log(min(c(sapply(res.ref,function(x){x$EndLambda}),recursive=T)))
    stopifnot(!is.na(lmin))
    lmax = log(max(c(sapply(res.ref,function(x){x$BeginLambda}),recursive=T)))
    lambdas = exp(seq(lmin,lmax,length.out=cv.gridsize))
    #build cv curve
    cv = foreach (lambda2=lambdas, .combine=rbind) %dopar% {#} %:% foreach (lambda1=lambdas, .combine=rbind) %dopar% {
      #compute coefficients on each model
      submat.all = mat[,.(groupname,ibin1,ibin2,ncounts)]
      submat.all[,value.ref:=flsaGetSolution(res.ref[[groupname[1]]], lambda1=0, lambda2=lambda2)[1,],by=groupname]
      #compute cv coefs for each group, given lambda
      submat.cv = foreach(cv=cvgroups, .combine=rbind) %do% {
        submat.cv = mat[,.(groupname,cv.group,ibin1,ibin2)]
        submat.cv[,value:=flsaGetSolution(res.cv[[cv.group[1]]][[groupname[1]]],
                                          lambda1=0, lambda2=lambda2)[1,],by=groupname]
        submat.cv=submat.cv[cv.group==cv]
      }
      #ggplot(submat.cv)+geom_raster(aes(ibin1,ibin2,fill=value))+facet_grid(groupname~cv.group)
      ggplot(submat.all)+geom_raster(aes(ibin1,ibin2,fill=value.ref))+facet_grid(~groupname)
      ret=merge(submat.cv, submat.all, all.x=T, by=c("groupname","ibin1","ibin2"), sort=F)
      ret=ret[,.(mse=mean((value-value.ref)^2)*ncounts[1]),by=c("groupname","ibin1","ibin2")]
      ret=ret[,.(mse=mean(mse)),by=groupname]
      ret[,.(groupname,lambda1=0,lambda2=lambda2,mse)]
    }
    ggplot(cv[lambda1==0&lambda2>7])+geom_line(aes(lambda2,mse,group=log(lambda1)))+facet_wrap(~groupname,scales="free")+scale_y_log10()
    ggplot(cv)+geom_raster(aes(log(lambda1),log(lambda2),fill=mse))+facet_wrap(~groupname)
    
  }
  counts[,c("bin1","bin2","ibin1","ibin2"):=list(NULL,NULL,NULL,NULL)]
  biases[,c("bin","ibin"):=list(NULL,NULL)]
  mat
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
group_datasets = function(cs, resolution, group=c("condition","replicate","enzyme","experiment"), ice=-1, verbose=T, ncores=1) {
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

#' Detect significant interactions wrt expected
#' 
#' @param cs CSnorm object
#' @param resolution,group see
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to
#'   identify the input matrices.
#' @param threshold significance threshold, between 0 and 1
#' @param ncores number of cores used for parallelization
#' @param use fused lasso detection (default FALSE)
#' @inheritParams call_interactions
#'   
#' @return the binned matrix with additional information relating to these 
#'   significant interactions
#' @export
#' 
#' @examples
detect_interactions = function(cs, resolution, group, threshold=0.95, ncores=1, binless=F){
  #get CSmat object
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (binless==T) {
    if (get_cs_interaction_idx(csm, type="binteractions", threshold=threshold, ref="expected", raise=F)>0) {
      stop("Refusing to overwrite this already detected interaction")
    }
    mat = csnorm_detect_binless(cs, resolution=resolution, group=group, ref="expected", threshold=threshold, ncores=ncores)
    csi=new("CSinter", mat=mat, type="binteractions", threshold=threshold, ref="expected")
  } else {
    if (get_cs_interaction_idx(csm, type="interactions", threshold=threshold, ref="expected", raise=F)>0) {
      stop("Refusing to overwrite this already detected interaction")
    }
    mat = csnorm_detect_binned(cs, resolution=resolution, group=group, ref="expected", threshold=threshold, ncores=ncores)
    csi=new("CSinter", mat=mat, type="interactions", threshold=threshold, ref="expected")
  }
  #store back
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
  mat=csnorm_detect_binned(cs, resolution=resolution, group=group, ref=ref, threshold=threshold, ncores=ncores)
  #add interaction to cs object
  csi=new("CSinter", mat=mat, type="differences", threshold=threshold, ref=ref)
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}
  
