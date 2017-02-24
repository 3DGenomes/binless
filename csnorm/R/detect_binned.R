#' @include csnorm.R
NULL

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

#' Binned detection of significant interactions wrt expected
#' 
#' @param cs CSnorm object
#' @param resolution,group see
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to
#'   identify the input matrices.
#' @param threshold significance threshold, between 0 and 1
#' @param ncores number of cores used for parallelization
#'   
#' @return the binned matrix with additional information relating to these 
#'   significant interactions
#' @export
#' 
#' @examples
detect_binned_interactions = function(cs, resolution, group, threshold=0.95, ncores=1){
  #get CSmat object
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csm, type="interactions", threshold=threshold, ref="expected", raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat = csnorm_detect_binned(cs, resolution=resolution, group=group, ref="expected", threshold=threshold, ncores=ncores)
  csi=new("CSinter", mat=mat, type="interactions", threshold=threshold, ref="expected")
  #store back
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}


#' Binned detection of significant differences with a reference
#' 
#' @inheritParams detect_binned_interactions
#'   
#' @return the binned matrix with additional information relating to these
#'   significant interactions
#' @export
#' 
#' @examples
detect_binned_differences = function(cs, resolution, group, ref, threshold=0.95, ncores=1){
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  if (get_cs_interaction_idx(csm, type="differences", threshold=threshold, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat = csnorm_detect_binned(cs, resolution=resolution, group=group, ref=ref, threshold=threshold, ncores=ncores)
  csi=new("CSinter", mat=mat, type="differences", threshold=threshold, ref=ref)
  #add interaction to cs object
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}
  
