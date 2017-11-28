#' @include binless.R
NULL

#' Prepare grouped signal matrix and settings
#' @keywords internal
prepare_signal_estimation = function(cs, csg, resolution, tol.val) {
  mat = binless:::get_signal_matrix(cs, resolution, groups=csg@names)
  #
  #add trail information
  stopifnot(all(mat[,.N,by=name]$N==mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  diag.rm = ceiling(cs@settings$dmin/resolution)
  #add other settings
  settings=list(metadata = get_signal_metadata(cs, csg@cts, resolution),
                nbins = csg@par$nbins,
                dispersion = csg@par$alpha,
                tol.val = tol.val,
                nperf = 50,
                min.patchsize = 4,
                min.l10FC = 0.5)
  cts=csg@cts[,.(name,bin1,bin2,count,lmu.nosig,weight,log_decay)]
  csi=new("CSbsig", mat=mat, cts=cts, settings=settings)
  return(csi)
}

#' Build grouped difference matrix using normalization data if available
#' @keywords internal
prepare_difference_estimation = function(cs, csg, resolution, ref, tol.val) {
  csi = binless:::prepare_signal_estimation(cs, csg, resolution, tol.val)
  names=csg@names
  mat = foreach(n=names[groupname!=ref,unique(groupname)],.combine=rbind) %do%
    merge(csi@mat[name==n],csi@mat[name==ref,.(bin1,bin2,phi1=phi)],all=T,by=c("bin1","bin2"))
  mat[,c("phi.ref","delta"):=list(phi1,(phi-phi1)/2)]
  mat[,c("phi","phi1"):=NULL]
  cts=csg@cts[name!=ref,.(name,bin1,bin2,count,lmu.nosig,weight,log_decay)]
  cts.ref=csg@cts[name==ref,.(name,bin1,bin2,count,lmu.nosig,weight,log_decay)]
  csi=new("CSbdiff", mat=mat, cts=cts, cts.ref=cts.ref,
          ref=as.character(ref), settings=csi@settings)
  return(csi)
}

#' Perform binless interaction detection using fused lasso
#'
#' @param cs 
#' @param ref 
#' @param resolution 
#' @param group 
#' @param ncores 
#' @param niter number of IRLS iterations, and BIC iterations within
#'   
#' @return 
#' @export
#' 
#' @examples
detect_binless_interactions = function(cs, resolution=cs@settings$base.res, group="all", ncores=1, tol.val=cs@settings$tol, verbose=T,
                                       fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA){
  if (verbose==T) cat("Binless interaction detection with resolution=",resolution," and group=",group,"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csg, type="CSbsig", raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  #
  ### prepare signal estimation
  if (verbose==T) cat("  Prepare for signal estimation\n")
  csi = binless:::prepare_signal_estimation(cs, csg, resolution, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  csigs = foreach(g=groupnames) %do% {
    csig = new("CSbsig", mat=csi@mat[name==g], cts=csi@cts[name==g], settings=csi@settings)
    csig@settings$last.beta=csi@mat[name==g,phi]
    csig@state = binless:::gfl_compute_initial_state(csig, diff=F)
    csig
  }
  registerDoParallel(cores=min(ncores,length(groupnames)))
  params = foreach(csig=csigs, .combine=rbind) %dopar% {
    binless:::fused_lasso(csig, positive=T, fixed=T, constrained=F, verbose=verbose,
                                fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                                fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at)
  }
  stopImplicitCluster()
  #display param info
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],"\n")
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  ### store interaction
  #store back
  csi@par=list(lambda1=params[,lambda1],lambda2=params[,lambda2],eCprime=params[,eCprime],name=params[,name])
  mat[,ncounts:=NULL]
  csi@mat=mat
  csg@interactions=append(csg@interactions,list(csi))
  cs@groups[[idx1]]=csg
  return(cs)
}


#' Binless detection of significant differences with a reference
#' 
#' @inheritParams detect_binless_interactions
#'   
#' @return the binned matrix with additional information relating to these
#'   significant interactions
#' @export
#' 
#' @examples
detect_binless_differences = function(cs, ref, resolution=cs@settings$base.res, group="all", ncores=1, tol.val=cs@settings$tol, verbose=T,
                                      fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=", group," and ref=",as.character(ref),"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="CSbdiff", ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  if (is.character(ref)) ref=csg@names[as.character(groupname)==ref,unique(groupname)]
  if (length(ref) == 0) stop(paste0("Invalid ref! For group by ",group,", acceptable refs are: ",paste(levels(ref),collapse=" / ")))
  if (verbose==T) cat("  Prepare for difference estimation\n")
  csi = binless:::prepare_difference_estimation(cs, csg, resolution, ref, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  csigs = foreach (g=groupnames) %do% {
    csig = new("CSbdiff", mat=csi@mat[name==g], cts=csi@cts[name==g], cts.ref=csi@cts.ref,
               ref=csi@ref, settings=csi@settings)
    csig@settings$last.beta=csi@mat[name==g,delta]
    csig@settings$last.phi.ref=csi@mat[name==g,phi.ref]
    csig@state = binless:::gfl_compute_initial_state(csig, diff=T)
    csig
  }
  registerDoParallel(cores=min(ncores,length(groupnames)))
  params = foreach(csig=csigs, .combine=rbind) %do% {
    binless:::fused_lasso(csig, positive=F, fixed=T, constrained=F, verbose=verbose,
                                fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                                fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at)
  }
  stopImplicitCluster()
  #display param info
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],"\n")
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  #store back
  csi@par=list(lambda1=params[,lambda1],lambda2=params[,lambda2],eCprime=params[,eCprime],name=params[,name])
  mat[,ncounts:=NULL]
  csi@mat=mat
  csg@interactions=append(csg@interactions,list(csi))
  cs@groups[[idx1]]=csg
  return(cs)
}
  
