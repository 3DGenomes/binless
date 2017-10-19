#' @include binless.R
NULL

#' Binned detection of significant interactions wrt expected
#' 
#' @param cs CSnorm object
#' @param resolution,group see
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to
#'   identify the input matrices.
#' @param ncores number of cores used for parallelization
#' @param threshold on the probability K/(1+K) where K is the Bayes factor
#' @param niter integer. Maximum number of IRLS iterations
#' @param tol numeric. Convergence tolerance for IRLS objective
#' @param verbose boolean.
#'   
#' @return the binned matrix with additional information relating to these 
#'   significant interactions
#' @export
#' 
#' @examples
detect_binned_interactions = function(cs, resolution, group, threshold=0.95, ncores=1,
                                      niter=100, tol=1e-3, verbose=T){
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csg, type="CSsig", threshold=threshold, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  cts = copy(csg@cts)
  #
  ### interaction detection
  if (verbose==T) cat("   Interaction detection\n")
  cts[,c("signal","mu.nosig"):=list(1,exp(lmu.nosig))]
  prior.sd=10
  for (i in 1:niter) {
    cts[,c("z","var","signal.old"):=list(count/(signal*mu.nosig)-1,
                                         (1/(signal*mu.nosig)+1/csg@par$alpha),signal)]
    cts[,phihat:=weighted.mean(z+log(signal), weight/var),by=c("name","bin1","bin2")]
    cts[,sigmasq:=1/sum(weight/var),by=c("name","bin1","bin2")]
    cts[,signal:=exp(phihat/(1+sigmasq/prior.sd^2))]
    if(cts[,all(abs(signal-signal.old)<tol)]) break
  }
  if (i==niter) cat("Warning: Maximum number of IRLS iterations reached for signal estimation!\n")
  #
  ### put in matrix form
  mat = cts[,.(K=exp(dnorm(phihat[1],mean=0,sd=sqrt(sigmasq[1]+prior.sd^2), log=T)-
                       dnorm(phihat[1],mean=0,sd=sqrt(sigmasq[1]), log=T)),
               signal=signal[1],
               signal.sd=sqrt(1/(1/sigmasq[1]+1/prior.sd^2))*signal[1],
               binned=sum(count)),keyby=c("name","bin1","bin2")]
  mat[,direction:=ifelse(signal>=1,"enriched","depleted")]
  mat[binned==0,c("signal","signal.sd","K","direction"):=list(0,NA,0,"depleted")]
  mat[,prob.gt.expected:=K/(1+K)]
  mat[,c("K","binned"):=list(NULL,NULL)]
  mat[,is.significant:=prob.gt.expected > threshold]
  #
  ### store back
  csi=new("CSsig", mat=mat, threshold=threshold)
  csg@interactions=append(csg@interactions,list(csi))
  cs@groups[[idx1]]=csg
  return(cs)
}


#' Binned detection of significant differences with a reference
#' 
#' @param ref character. The name of a reference group for difference detection.
#' @inheritParams detect_binned_interactions
#'   
#' @return the binned matrix with additional information relating to these
#'   significant interactions
#' @export
#' 
#' @examples
detect_binned_differences = function(cs, resolution, group, ref, threshold=0.95, ncores=1,
                                     niter=100, tol=1e-3, verbose=T){
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  #
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csg, type="CSdiff", threshold=threshold, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  # check if reference is ok
  cts = csg@cts
  if (!(ref %in% cts[,unique(name)]))
    stop("Reference group name not found! Valid names are: ",deparse(as.character(cts[,unique(name)])))
  if (cts[name!=ref,.N]==0)
    stop("There is no other group than ",ref, ", cannot compute differences!")
  #
  if (verbose==T) cat("   Difference detection\n")
  #replicate reference counts for each case
  cts[,c("mu.nosig"):=list(exp(lmu.nosig))]
  ctsref = foreach(n=cts[name!=ref,unique(name)],.combine=rbind) %do%
    cts[name==ref,.(name=n,bin1,bin2,count,mu,weight,mu.nosig)]
  cts=cts[name!=ref]
  #IRLS iteration
  mat=cts[,.(phi.ref=0,delta=0,diffsig=1),by=c("name","bin1","bin2")]
  prior.sd=10
  for (i in 1:niter) {
    if (i>1) {
      mat[,phi.ref:=(phihat.ref/sigmasq.ref + (phihat-delta)/sigmasq)/(1/sigmasq.ref+1/sigmasq)]
      ctsref=ctsref[,.(name,bin1,bin2,count,mu.nosig,weight)]
      cts=cts[,.(name,bin1,bin2,count,mu.nosig,weight)]
      mat=mat[,.(name,bin1,bin2,phi.ref,delta,diffsig,diffsig.old)]
    }
    ctsref=mat[ctsref]
    ctsref[,c("z","var"):=list(count/(exp(phi.ref)*mu.nosig)-1,(1/(exp(phi.ref)*mu.nosig)+1/csg@par$alpha))]
    mat=mat[ctsref[,.(phihat.ref=weighted.mean(z+phi.ref, weight/var),
                      sigmasq.ref=1/sum(weight/var)),keyby=c("name","bin1","bin2")]]
    #
    cts=mat[cts]
    cts[,c("z","var"):=list(count/(exp(phi.ref+delta)*mu.nosig)-1,
                            (1/(exp(phi.ref+delta)*mu.nosig)+1/csg@par$alpha))]
    mat=mat[cts[,.(phihat=weighted.mean(z+phi.ref+delta, weight/var),
                   sigmasq=1/sum(weight/var)),by=c("name","bin1","bin2")]]
    mat[,deltahat:=phihat-phihat.ref]
    mat[,delta:=deltahat/(1+(sigmasq.ref+sigmasq)/prior.sd^2)]
    mat[,diffsig.old:=diffsig]
    mat[,diffsig:=exp(delta)]
    if(mat[,all(abs(diffsig-diffsig.old)<tol)]) break
  }
  if (i==niter) message("Warning: Maximum number of IRLS iterations reached for signal estimation!\n")
  #
  # remove extra columns and compute Bayes factor
  mat=mat[,.(name,bin1,bin2,difference=diffsig,
             delta.sd=(sigmasq.ref+sigmasq)/(1+(sigmasq.ref+sigmasq)/prior.sd^2),
             K=exp(dnorm(deltahat,mean=0,sd=sqrt(sigmasq+sigmasq.ref+prior.sd^2),log=T)-
                     dnorm(deltahat,mean=0,sd=sqrt(sigmasq+sigmasq.ref),log=T)))]
  mat[,difference.sd:=delta.sd*difference]
  mat[,delta.sd:=NULL]
  # treat limiting cases
  mat=cts[,.(is.zero=sum(count)==0),by=c("name","bin1","bin2")][mat]
  mat=ctsref[,.(is.zero.ref=sum(count)==0),by=c("name","bin1","bin2")][mat]
  mat[is.zero==T,c("difference","difference.sd","K"):=list(0,NA,0)]
  mat[is.zero.ref==T,c("difference","difference.sd","K"):=list(Inf,NA,Inf)]
  mat[is.zero.ref==T&is.zero==T,c("difference","difference.sd","K"):=list(NA,NA,NA)]
  mat[,c("is.zero","is.zero.ref"):=NULL]
  # additional columns
  colname=paste0("prob.gt.",ref)
  mat[,c(colname):=K/(1+K)]
  mat[,is.significant:=get(colname) > threshold]
  mat[,c("direction","K"):=list(ifelse(difference>=1,"enriched","depleted"),NULL)]
  #
  ### store back
  csi=new("CSdiff", mat=mat, threshold=threshold, ref=as.character(ref))
  csg@interactions=append(csg@interactions,list(csi))
  cs@groups[[idx1]]=csg
  return(cs)
}
  
