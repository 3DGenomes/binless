#' @include csnorm.R
NULL

#' Convert a sparse counts data.table to a dense one by adding rows with zero counts
#'
#' @param counts,biases data.tables as returned by \code{\link{prepare_for_sparse_cs_norm}}
#' @param biases2 data.table of biases for id2 column of counts. If NULL (default), use that of biases
#'
#' @return a counts data.table with zeros filled according to cut sites provided in biases (and biases2 if available)
#' @keywords internal
#' @export
#' @section Warning:
#' Memory-intensive
#' @examples
fill_zeros = function(counts,biases,biases2=NULL,circularize=-1L) {
  if (is.null(biases2)) biases2=biases
  runname=biases[,unique(name)]
  if (length(runname)>1) {
    foreach (i=runname, .combine=rbind) %do%
      fill_zeros(counts[name==i],biases[name==i],biases2[name==i],circularize=circularize)
  } else {
    if (biases[,.N]==0 | biases2[,.N]==0) return(data.table())
    newcounts=CJ(biases[,paste(id,pos)],biases2[,paste(id,pos)])
    newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
    newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
    newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=
                list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
    newcounts=newcounts[pos1<pos2]
    setkey(newcounts, id1, id2, pos1, pos2)
    setkey(counts, id1, id2, pos1, pos2)
    newcounts=counts[newcounts]
    newcounts[is.na(contact.close),contact.close:=0]
    newcounts[is.na(contact.far),contact.far:=0]
    newcounts[is.na(contact.up),contact.up:=0]
    newcounts[is.na(contact.down),contact.down:=0]
    newcounts[is.na(name),name:=runname]
    if (circularize>0) {
      newcounts[,distance:=pmin(abs(pos2-pos1), circularize+1-abs(pos2-pos1))]
    } else {
      newcounts[,distance:=abs(pos2-pos1)]
    }
    newcounts
  }
}


#' Single-cpu fitting
#' @keywords internal
#' 
csnorm_fit = function(biases, counts, design, dmin, dmax, bf_per_kb=1, bf_per_decade=5, iter=10000,
                      verbose=T, init=0, weight=array(1,dim=design[,.N]), ...) {
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=design[,.N], Biases=design[,uniqueN(genomic)], Decays=design[,uniqueN(decay)],
               XB=as.array(design[,genomic]), XD=as.array(design[,decay]),
               Krow=Krow, SD=biases[,.N], bbegin=bbegin,
               cutsitesD=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cbegin=cbegin,
               cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight))
  op=optimizing(stanmodels$fit, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  return(op)
}

#' Predict expected values for each count given optimized model parameters
#' 
#' @param cs CSnorm object
#' @param counts
#' @param verbose
#'   
#' @section Warning: Do not call this function to compute a binned matrix. It
#'   should (almost) never be used and might be removed in the future.
#'   
#' @return a stan optimization output with the predictions.
#' @keywords internal
#' @export
#' 
#' @examples
csnorm_predict_all = function(cs, counts, verbose=T) {
  biases=cs@biases
  par=cs@par
  dmin=cs@settings$dmin
  dmax=cs@settings$dmax
  Kdiag=round((log10(dmax)-log10(dmin))*cs@settings$bf_per_decade)
  ncounts=counts[,.N]
  if (counts[,.N]==1) {
    cbegin=c(1,2)
  } else {
    cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  }
  design=cs@design
  data = list( Dsets=design[,.N], Decays=design[,uniqueN(decay)], XD=as.array(design[,decay]),
               Kdiag=Kdiag, SD=biases[,.N], cutsitesD=biases[,pos], dmin=dmin, dmax=dmax,
               N=counts[,.N], cbegin=cbegin, cidx=t(data.matrix(counts[,.(id1,id2)])), dist=as.array(counts[,distance]),
               eC=par$eC, log_nu=par$log_nu, log_delta=par$log_delta,
               beta_diag_centered=par$beta_diag_centered)
  capture.output(pred<-as.data.table(optimizing(stanmodels$predict_all, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)$par))
  pred=cbind(counts,pred)
  return(pred)
}

#' Fill in a large matrix with predicted values
#' @keywords internal
#' @export
#'
csnorm_predict_all_parallel = function(cs, counts, verbose=T, ncores=1) {
  nchunks=min(10*ncores,counts[,.N,by=name][,N]) #ensure at least 1 of each per parallel prediction
  counts[,chunk:=.I]
  counts[,chunk:=chunk-min(chunk),by=name]
  counts[,chunk:=as.integer(chunk/((max(chunk)+1)/nchunks)),by=name]
  registerDoParallel(cores = ncores)
  counts = foreach (i=0:(nchunks-1), .combine=rbind) %dopar%
    csnorm_predict_all(cs,counts[chunk==i],verbose=verbose)
  counts[,chunk:=NULL]
  counts
}

#' Run exact model on a single cpu
#' @inheritParams run_exact
#' @keywords internal
#' @export
#' 
run_serial = function(cs, bf_per_kb=1, bf_per_decade=5, lambda=1, iter=100000, subsampling.pc=100) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, lambda=lambda, iter=iter))
  #fill counts matrix
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize)
  #report min/max distance
  dmin=0.99
  if (cs@settings$circularize>0) {
    dmax=cs@settings$circularize/2+0.01
  } else {
    dmax=cs@biases[,max(pos)-min(pos)]+0.01
  }
  cs@settings$dmin=dmin
  cs@settings$dmax=dmax
  setkey(cs@biases,name,id,pos)
  setkey(cs@counts,name,id1,pos1,id2,pos2)
  #initial guess
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases, design=cs@design,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, lambda=lambda, verbose=T, iter=iter)))
  cs@diagnostics=list(out.init=init.output, runtime.init=init.a[1]+init.a[4])
  #main optimization, subsampled
  counts.sub=cs@counts[sample(.N,round(subsampling.pc/100*.N))]
  setkeyv(counts.sub,key(cs@counts))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = counts.sub, design=cs@design, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=iter, verbose = T,
    init=init.op, weight=counts.sub[,.N,by=name]$N/cs@counts[,.N,by=name]$N)))
  cs@diagnostics=list(out=output, runtime=a[1]+a[4])
  #report statistics
  op$par$init=init.op
  op$par$value=op$value
  if (subsampling.pc<100) op$par$counts.sub=counts.sub
  cs@par=op$par
  cs
}

#' Cut-site normalization (exact model)
#' 
#' Will run the exact model of normalization (on one cpu for each lambda 
#' provided) and returns the most likely model and predicted quantities. Useful
#' for comparison purposes. If you don't know what to use, try 
#' \code{\link{run_simplified}}.
#' 
#' @inheritParams run_simplified
#' @param subsampling.pc numeric. Percentage of the data used to do the calculations (default 100).
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
run_exact = function(cs, bf_per_kb=1, bf_per_decade=5, lambdas=c(0.1,1,10), ncores=1, iter=100000, subsampling.pc=100) {
  cs@binned=list() #erase old binned datasets if available
  registerDoParallel(cores=ncores)
  cs = foreach (lambda=lambdas, .combine=function(x,y){if (x@par$value<y@par$value){return(y)}else{return(x)}}) %dopar%
    run_serial(cs, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, lambda=lambda, iter=iter, subsampling.pc=subsampling.pc)
  return(cs)
}

#' Verify model fit by computing posterior predictive quantities
#'
#' @param cs CSnorm object, normalized.
#' @param genomic.groups How many groups for nu and delta
#' @param decay.groups How many groups for diagonal decay
#' @param npoints Number of points per group to take (not guaranteed)
#'
#' @return
#' @export
#'
#' @examples
check_fit = function(cs, genomic.groups=5, decay.groups=5, npoints=10) {
  if (length(cs@par)==0) stop("Must normalize the datasets first")
  #build bins
  dbins=c(0,10**seq(3,log10(cs@settings$dmax),length.out=decay.groups))
  gbins=cut2(c(cs@par$log_nu,cs@par$log_delta),g=genomic.groups,onlycuts=T)
  #build counts matrix 
  biases=copy(cs@biases[,.(name,id,pos)])
  biases[,c("log_nu","log_delta"):=list(cs@par$log_nu,cs@par$log_delta)]
  biases[,nubin:=cut(log_nu, gbins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  biases[,deltabin:=cut(log_delta, gbins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  biases[,c("log_nu","log_delta"):=list(NULL,NULL)]
  biases=biases[,.SD[sample(.N,min(.N,npoints))],by=c("name","nubin","deltabin")]
  counts=cs@counts[id1%in%biases[,id]&id2%in%biases[,id]]
  counts=fill_zeros(counts,biases,circularize=cs@settings$circularize)
  #filter by distance
  counts[,dbin:=cut(distance, dbins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  counts=counts[,.SD[sample(.N,min(.N,npoints))],by=c("name","dbin")]
  setkey(counts, name, id1, id2, pos1, pos2)
  #predict values
  counts=csnorm_predict_all(cs, counts, verbose=F)
  counts=rbind(counts[,.(name,id1,id2,dbin,count=contact.close,mean=exp(log_mean_cclose))],
               counts[,.(name,id1,id2,dbin,count=contact.far,mean=exp(log_mean_cfar))],
               counts[,.(name,id1,id2,dbin,count=contact.up,mean=exp(log_mean_cup))],
               counts[,.(name,id1,id2,dbin,count=contact.down,mean=exp(log_mean_cdown))])
  biases[,pos:=NULL]
  counts=merge(counts,biases,by.x=c("name","id1"),by.y=c("name","id"))
  counts=merge(counts,biases,by.x=c("name","id2"),by.y=c("name","id"),suffixes=c("1","2"))
  #compute p-values
  counts[count>=mean,pval:=pnbinom(count,size=cs@par$alpha,mu=mean,lower.tail=F)]
  counts[count<mean,pval:=pnbinom(count,size=cs@par$alpha,mu=mean,lower.tail=T)]
  counts[,sd:=sqrt(mean+mean**2/cs@par$alpha)]
  #graph p-values
  p.all=ggplot(counts)+geom_histogram(aes(pval))+facet_wrap(~name)+xlab("model p-value")+ylab("frequency")
  p.decay=ggplot(counts)+geom_jitter(aes(1,pval))+facet_grid(name~dbin)+ylab("model p-value")+xlab("diagonal decay bin")
  p.nu=ggplot(counts)+geom_jitter(aes(1,pval))+facet_grid(name~nubin1)+ylab("model p-value")+xlab("nu bin")
  p.delta=ggplot(counts)+geom_jitter(aes(1,pval))+facet_grid(name~deltabin1)+ylab("model p-value")+xlab("delta bin")
  return(list(all=p.all,decay=p.decay,nu=p.nu,delta=p.delta,counts=counts))
}
