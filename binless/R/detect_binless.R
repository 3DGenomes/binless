#' @include binless.R
NULL

#' compute BIC for a given value of lambda2, optimizing lambda1 and eCprime (performance iteration, persistent state)
#' @keywords internal
gfl_BIC = function(csig, lambda2, constrained=T, positive=T, fixed=F, fix.lambda1=F, fix.lambda1.at=NA) {
  stopifnot(fix.lambda1==F || fix.lambda1.at>=0) #otherwise eCprime cannot be determined
  stopifnot(fixed==T || positive==T) #the case positive==F && fixed==F is not implemented
  #
  ctsg=csig@cts
  nbins=csig@settings$nbins
  dispersion=csig@settings$dispersion
  metadata=csig@settings$metadata
  tol.val=csig@settings$tol.val
  state=csig@state
  stopifnot(length(state)>0) #warm start only
  if (class(csig)=="CSbdiff") ctsg.ref=csig@cts.ref else ctsg.ref=NULL
  nperf=csig@settings$nperf
  if (is.null(ctsg.ref)) {
    perf.c = binless:::wgfl_signal_BIC(ctsg, dispersion, nperf, nbins, state$GFLState, lambda2,
                                      tol.val, metadata, state$beta, constrained, fixed, positive, fix.lambda1, fix.lambda1.at)
  } else {
    perf.c = binless:::wgfl_diff_BIC(ctsg, ctsg.ref, dispersion, nperf, nbins, state$GFLState, lambda2,
                                      tol.val, metadata, state$phi.ref, state$beta, constrained, fix.lambda1, fix.lambda1.at)
  }
  return(perf.c)
}

#' cross-validate lambda2 assuming it is very rugged
#' @keywords internal
optimize_lambda2 = function(csig, n.SD=1, constrained=T, positive=T, fixed=F, fix.lambda1=F, fix.lambda1.at=NA, ncores=1) {
  obj = function(x) {
    csig@state <<- binless:::gfl_BIC(csig, lambda2=10^(x), constrained=constrained, positive=positive, fixed=fixed, fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
    #cat("optimize_lambda2: eval at lambda2= ",csig@state$lambda2, " lambda1= ",csig@state$lambda1,
    #    " eCprime= ",csig@state$eCprime," BIC= ",csig@state$BIC, " dof= ",csig@state$dof,"\n")
    return(csig@state$BIC)
  }
  #
  # dt3 = foreach (lam=10^seq(-1,1,length.out=20),.combine=rbind) %:% foreach(csig=csigs,.combine=rbind) %dopar% {
  #   dset = csig@cts[,name[1]]
  #   state = binless:::gfl_BIC(csig, lambda2=lam, constrained=constrained, positive=positive, fixed=fixed, fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
  #   csig@state = state
  #   mat=as.data.table(state$mat)
  #   state$dset=dset
  #   info=as.data.table(state[c("lambda2","lambda1","eCprime","dof","BIC","BIC.sd","converged","dset")])
  #   cbind(mat,info)
  # }
  #ggplot(dt3[,.SD[1],by=c("lambda2","dset")])+geom_point(aes(lambda2,BIC,colour=converged,group=dset))+facet_grid(dset~., scales="free_y")#+scale_y_log10()
  #ggplot(dt3[,.SD[1],by=c("lambda2","dset")])+geom_line(aes(lambda2,BIC,colour=converged,group=dset))+geom_errorbar(aes(lambda2,ymin=BIC-BIC.sd,ymax=BIC+BIC.sd,colour=converged))+facet_grid(dset~., scales="free_y")#+scale_y_log10()
  #ggplot(dt3[,.SD[BIC==min(BIC)],by=c("dset")])+geom_raster(aes(bin1,bin2,fill=pmin(phihat,3)))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~dset)+coord_fixed()+scale_fill_gradient2()
  #ggplot(dt3[,.SD[abs(lambda2-1)==min(abs(lambda2-1))],by=c("dset")])+geom_raster(aes(bin1,bin2,fill=pmin(phihat,3)))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~dset)+coord_fixed()+scale_fill_gradient2()
  # ggplot(dt3[dset==levels(dset)[1]])+geom_raster(aes(bin1,bin2,fill=pmin(phihat,3)))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~factor(lambda2))+coord_fixed()+scale_fill_gradient2()
   
  # scores=rbind(dt3all[resolution=="5k",.SD[1],by=c("lambda2","dset","resolution")],
  #              dt3.5k.new[,.SD[1],by=c("lambda2","dset","resolution")])
  # ggplot(scores)+geom_line(aes(lambda2,BIC,colour=correct))+
  #   geom_hline(aes(yintercept=BIC+BIC.sd),data=scores[,.SD[BIC==min(BIC)],by=c("dset","resolution")])+
  #   geom_errorbar(aes(lambda2,ymin=BIC-BIC.sd,ymax=BIC+BIC.sd,colour=converged))+facet_wrap(dset~resolution, scales="free_y")#+scale_y_log10()
   
  # scores=dt3all[,.SD[1],by=c("lambda2","dset","resolution")]
  # scores[,upper:=BIC+BIC.sd]
  # ggplot(scores)+geom_line(aes(lambda2,BIC,colour=correct))+
  #   geom_hline(aes(yintercept=upper),data=scores[,.SD[BIC==min(BIC)],by=c("dset","resolution")])+
  #   geom_errorbar(aes(lambda2,ymin=BIC-BIC.sd,ymax=BIC+BIC.sd,colour=converged))+facet_wrap(~resolution+dset, scales="free")#+scale_y_log10()
  
  #ggplot(dcast(dt3,lambda2+bin1+bin2~ori,value.var = "phi"))+geom_raster(aes(bin1,bin2,fill=new))+geom_raster(aes(bin2,bin1,fill=old))+facet_wrap(~lambda2)+coord_fixed()+scale_fill_gradient2()
  #dt.fix = foreach (lam=10^(seq(log10(12.8),log10(13.3),length.out=20)),.combine=rbind) %dopar% {
  #  csig@state <<- binless:::gfl_BIC(csig, lambda2=lam, constrained=constrained, positive=positive, fixed=fixed)
  #  as.data.table(csig@state[c("lambda2","lambda1","eCprime","dof","BIC","BIC.sd","converged")])
  #}
  #dt=rbindlist(list(free=dt.free,fix=dt.fix), use=T, idcol="ori")
  #ggplot(dt)+geom_line(aes(lambda2,BIC,colour=ori))#+geom_point(data=r,aes(10^x,y))#+xlim(1,3)+ylim(1600,2000)
  ### optimization in four stages
  #first, find rough minimum between 1 and 100
  minlambda=1
  maxlambda=100
  op<-optimize(obj, c(log10(minlambda),log10(maxlambda)), tol=0.1)
  lambda2=10^op$minimum
  #second, repeated gridding at finer scales
  range=10
  npoints=20
  reduction=3
  nrounds=3
  dt=data.table()
  for (i in 1:nrounds) {
    minlambda = max(minlambda,lambda2 - range/2)
    maxlambda = min(maxlambda,lambda2 + range/2)
    #cat("round ",i,": min=",minlambda," < lambda2=",lambda2, " < max=",maxlambda," range=",range,"\n")
    dt = rbind(dt,foreach (lam=log10(seq(minlambda, maxlambda, l=npoints)),
                   .combine=rbind) %do% data.table(x=lam,y=obj(lam),i=paste(i)))
    lambda2 = dt[y==min(y),10^x]
    range=range/reduction
  }
  #ggplot()+geom_point(aes(10^x,y,colour=i),data=dt)
  #third, minimize fully around that minimum
  minlambda = max(minlambda,lambda2 - range/2)
  maxlambda = min(maxlambda,lambda2 + range/2)
  #cat("final round : min=",minlambda," < lambda2=",lambda2, " < max=",maxlambda," range=",range,"\n")
  op<-optimize(obj, log10(c(minlambda,maxlambda)))
  lambda2=10^op$minimum
  if (n.SD>0) {
    #finally, find minimum + SD
    minlambda=lambda2
    maxlambda=100
    optBIC=csig@state$BIC+n.SD*csig@state$BIC.sd
    obj2 = function(x) {
      a=obj(x)
      return(a+2*abs(optBIC-a))
    }
    op<-optimize(obj2, c(log10(minlambda),log10(maxlambda)))
    lambda2=10^op$minimum
  }
  #finish
  if (lambda2==maxlambda) cat("   Warning: lambda2 hit upper boundary.\n")
  if (lambda2==minlambda) cat("   Warning: lambda2 hit lower boundary.\n")
  obj(log10(lambda2))
  retvals = as.list(csig@state)[c("lambda2","lambda1","eCprime","BIC","BIC.sd","dof")]
  if (fixed==T && abs(retvals$eCprime)>csig@settings$tol.val) cat("Warning: fixed = T but eCprime != 0\n") #only when fix.lambda1==F
  csig@par=modifyList(csig@par,retvals)
  return(csig)
}

#' cross-validate lambda2 assuming it is smooth
#' @keywords internal
optimize_lambda2_smooth = function(csig, n.SD=1, constrained=T, positive=T, fixed=F, fix.lambda1=F, fix.lambda1.at=NA) {
  obj = function(x) {
    csig@state <<- binless:::gfl_BIC(csig, lambda2=10^(x), constrained=constrained, positive=positive, fixed=fixed, fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
    #cat("optimize_lambda2: eval at lambda2= ",csig@state$lambda2, " lambda1= ",csig@state$lambda1,
    #    " eCprime= ",csig@state$eCprime," BIC= ",csig@state$BIC, " dof= ",csig@state$dof,"\n")
    return(csig@state$BIC)
  }
  #first, find rough minimum between 1 and 100
  minlambda=1
  maxlambda=100
  op<-optimize(obj, c(log10(minlambda),log10(maxlambda)), tol=0.1)
  lambda2=10^op$minimum
  if (n.SD>0) {
    #finally, find minimum + SD
    minlambda=lambda2
    maxlambda=100
    optBIC=csig@state$BIC+n.SD*csig@state$BIC.sd
    obj2 = function(x) {
      a=obj(x)
      return(a+2*abs(optBIC-a))
    }
    op<-optimize(obj2, c(log10(minlambda),log10(maxlambda)), tol=0.1)
    lambda2=10^op$minimum
  }
  #finish
  if (lambda2==maxlambda) cat("   Warning: lambda2 hit upper boundary.\n")
  if (lambda2==minlambda) cat("   Warning: lambda2 hit lower boundary.\n")
  obj(log10(lambda2))
  retvals = as.list(csig@state)[c("lambda2","lambda1","eCprime","BIC","BIC.sd","dof")]
  if (fixed==T && abs(retvals$eCprime)>csig@settings$tol.val) cat("Warning: fixed = T but eCprime != 0\n") #only when fix.lambda1==F
  csig@par=modifyList(csig@par,retvals)
  return(csig)
}

evaluate_at_lambda2 = function(csig, lambda2, constrained=T, positive=T, fixed=F, fix.lambda1=F, fix.lambda1.at=NA) {
  csig@state = binless:::gfl_BIC(csig, lambda2=lambda2, constrained=constrained, positive=positive, fixed=fixed,
                              fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
  retvals = as.list(csig@state)[c("lambda2","lambda1","eCprime","BIC","BIC.sd","dof")]
  if (fixed==T && abs(retvals$eCprime)>csig@settings$tol.val) cat("Warning: fixed = T but eCprime != 0\n") #only when fix.lambda1==F
  csig@par=modifyList(csig@par,retvals)
  return(csig)
}

#' run fused lasso on one dataset contained in matg, fusing 'valuehat' into
#' 'value'
#' 
#' @param csig the CSbsig or CSbdiff object to optimize
#' @param positive boolean. Constrain eCprime in order to force 'value' to be
#'   positive?
#' @param fixed boolean. Set eCprime=0 throughout ?
#' @param constrained boolean. Constrain lambda1 so that any diagonal that contains
#'   only one big patch be forced to have 'value'=0 ?
#' @param verbose boolean (default TRUE).
#'   
#'   finds optimal lambda1, lambda2 and eC using BIC.
#' @keywords internal
fused_lasso = function(csig, positive, fixed, constrained, verbose=T, fix.lambda1=F, fix.lambda1.at=0.1, fix.lambda2=F, fix.lambda2.at=NA) {
  if (fix.lambda2==F) {
    n.SD=ifelse(fixed==T,1,0)
    csig = binless:::optimize_lambda2_smooth(csig, n.SD=n.SD, constrained=constrained, positive=positive, fixed=fixed,
                                   fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
  } else {
    stopifnot(fix.lambda2.at>0)
    csig = binless:::evaluate_at_lambda2(csig, fix.lambda2.at, constrained=constrained, positive=positive, fixed=fixed,
                                     fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
  }
  csig@par$name=csig@cts[,name[1]]
  matg = as.data.table(csig@state$mat)
  setkey(matg,bin1,bin2)
  matg[,name:=csig@par$name]
  params = as.data.table(csig@par)
  params[,c("state","mat"):=list(list(csig@state),list(matg))]
  return(params)
}

#' Build (grouped) signal matrix using normalization data if available
#' @keywords internal
get_signal_matrix = function(cs, resolution=cs@settings$base.res, groups=cs@experiments[,.(name,groupname=name)]) {
  #build signal matrix
  if (cs@par$signal[,.N]==0) {
    sbins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
    signal.bins=unique(cut(c(sbins,head(sbins,n=-1)+resolution/2), sbins,
                           ordered_result=T, right=F, include.lowest=T,dig.lab=12))
    mat=CJ(name=groups[,groupname],bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
    mat[,phi:=0]
  } else {
    #report phi values for each group
    if (resolution!=cs@settings$base.res) {
      refmat=groups[cs@par$signal[,.(name,bin1,bin2,phi)]][
        !is.na(groupname),.(name=groupname,refbin1=bin1,refbin2=bin2,phi)][
          ,.(phi=mean(phi)),by=c("name","refbin1","refbin2")]
      #merge signal to new binning
      sbins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
      pos=head(sbins,n=-1)+resolution/2
      bins=unique(data.table(refbin=cut(pos, cs@settings$sbins,
                                        ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                             bin=cut(pos, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)))
      stopifnot(bins[,.N]==length(sbins)-1)
      mat=CJ(name=groups[,groupname],bin1=bins[,unique(bin)],bin2=bins[,unique(bin)],
             sorted=F,unique=F)[bin2>=bin1]
      mat=merge(mat,bins,by.x="bin1",by.y="bin",all.y=T)
      mat=merge(mat,bins,by.x="bin2",by.y="bin",all.y=T, suffixes=c("1","2"))
      mat=merge(mat,refmat,by=c("name","refbin1","refbin2"),all.x=T)
      mat[is.na(phi),phi:=0]
      mat=mat[,.(name,bin1,bin2,phi)]
    } else {
      mat=groups[cs@par$signal[,.(name,bin1,bin2,phi)]][!is.na(groupname),.(name=groupname,bin1,bin2,phi)]
    }
    mat=mat[,.(phi=mean(phi)),keyby=c("name","bin1","bin2")]
  }
  #ggplot(cs@par$signal)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
  return(mat)
}

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

#' build initial state from phi / delta
#' 
#' @keywords internal
gfl_compute_initial_state = function(csig, diff=F) {
  matg = csig@mat
  if (diff==F) {
    state = list(beta=matg[,phi], lambda1=0.05, eCprime=0, GFLState = list())
  } else {
    state = list(phi.ref=matg[,phi.ref], beta=matg[,delta], GFLState = list())
  }
  return(state)
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
detect_binless_interactions = function(cs, resolution, group, ncores=1, tol.val=cs@settings$tol, verbose=T,
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
  mat[,c("value","valuehat"):=list(phi,phihat)]
  #
  #if (verbose==T) cat(" Detect patches\n")
  #mat = foreach(g=groupnames, .combine=rbind) %do% {
  #  matg = mat[name==g]
  #  #cl = binless:::build_patch_graph_components(csi@settings$nbins, matg, csi@settings$tol.val)
  #  #matg[,c("patchno","value"):=list(factor(cl$membership),NULL)]
  #  matg[,value:=phi]
  #  matg = binless:::detect_binless_patches(matg, csi@settings)
  #  matg[,value:=NULL]
  #  matg
  #}
  #
  ### store interaction
  #store back
  csi@par=list(lambda1=params[,lambda1],lambda2=params[,lambda2],name=params[,name])
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
detect_binless_differences = function(cs, resolution, group, ref, ncores=1, tol.val=cs@settings$tol, verbose=T,
                                      fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=", group," and ref=",as.character(ref),"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="CSbdiff", ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  if (is.character(ref)) ref=csg@names[as.character(groupname)==ref,unique(groupname)]
  if (verbose==T) cat("  Prepare for difference estimation\n")
  csi = binless:::prepare_difference_estimation(cs, csg, resolution, ref, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  csigs = foreach (g=groupnames) %do% {
    csig = new("CSbdiff", mat=csi@mat[name==g], cts=csi@cts[name==g], cts.ref=csi@cts.ref,
               ref=csi@ref, settings=csi@settings)
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
  mat[,c("value","valuehat"):=list(delta,deltahat)]
  #
  #if (verbose==T) cat(" Detect patches\n")
  #mat = foreach(g=groupnames, .combine=rbind) %do% {
  #  matg = mat[name==g]
  #  #cl = binless:::build_patch_graph_components(csi@settings$nbins, matg, csi@settings$tol.val)
  #  #matg[,c("patchno","value"):=list(factor(cl$membership),NULL)]
  #  matg[,value:=delta]
  #  matg = binless:::detect_binless_patches(matg, csi@settings)
  #  matg[,value:=NULL]
  #  matg
  #}
  #store back
  csi@par=list(lambda1=params[,lambda1],lambda2=params[,lambda2],name=params[,name])
  csi@mat=mat
  csg@interactions=append(csg@interactions,list(csi))
  cs@groups[[idx1]]=csg
  return(cs)
}
  
#' make plot of binless matrix, with minima/maxima highlighted
#'
#' @param mat the binless matrix
#' @param minima whether to display minima or not (say TRUE for difference matrix) 
#'
#' @return
#' @export
#'
#' @examples
plot_binless_matrix = function(mat, minima=F, scale=T, value.name="value") {
  if (length(setdiff(c("begin1","begin2","end1","end2"),names(mat)))>0) {
    mat[,c("begin1","begin2","end1","end2"):=list(unclass(bin1),unclass(bin2),unclass(bin1)+1,unclass(bin2)+1)]
  }
  resolution=mat[bin1==bin1[1]&name==name[1],begin2[2]-begin2[1]]
  a=mat[is.maximum==T]
  a=a[,.SD[,.(begin1=c(begin1,begin1,end1,end1)-resolution/2, begin2=c(begin2,end2,begin2,end2)-resolution/2,
              patchno, get(value.name))][chull(begin1,begin2)], by=c("patchno","name")]
  if (minima==T) {
    b=mat[is.minimum==T]
    b=b[,.SD[,.(begin1=c(begin1,begin1,end1,end1)-resolution/2, begin2=c(begin2,end2,begin2,end2)-resolution/2,
                patchno, get(value.name))][chull(begin1,begin2)], by=c("patchno","name")]
    p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=(get(value.name))))+
      geom_raster(aes(begin2,begin1,fill=(get(value.name))))+
      scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value = "white")+facet_wrap(~name)+labs(fill=value.name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="blue",fill=NA,data=b)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="red",fill=NA,data=a)
  } else {
    p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=get(value.name)))+
      geom_raster(aes(begin2,begin1,fill=get(value.name)))+
      scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value = "white")+facet_wrap(~name)+labs(fill=value.name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="black",fill=NA,data=a)
  }
  if (scale!=T) p=p+guides(fill=F)
  print(p)
}