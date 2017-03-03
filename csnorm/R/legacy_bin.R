#' @include csnorm.R
NULL

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
