#' @include csnorm.R
NULL

#' dispatch for fused lasso perf iteration: warm/cold signal/difference
#' @keywords internal
gfl_perf_iteration = function(csig, lambda1, lambda2, eCprime) {
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
    perf.c = csnorm:::wgfl_signal_perf_warm(ctsg, dispersion, nperf, nbins, state$GFLState, lambda2, tol.val/20, metadata, state$beta)
  } else {
    stopifnot(eCprime==0)
    perf.c = csnorm:::wgfl_diff_perf_warm(ctsg, ctsg.ref, dispersion, nperf, nbins, state$GFLState, lambda2, tol.val/20, metadata,
                                            state$phi.ref, state$beta)
  }
  return(perf.c)
}

#' compute sparse fused lasso matrix for a given value of lambda1, lambda2 and eCprime (performance iteration)
#' @keywords internal
gfl_get_matrix = function(csig, lambda1, lambda2, eCprime) {
  if (csig@state$lambda1==lambda1 & csig@state$lambda2==lambda2 & (class(csig)=="CSbdiff" | csig@state$eCprime==eCprime)) {
    perf.c = csig@state
  } else {
    if (eCprime != 0 | lambda1>0) cat("WARNING: returned matrix will be unthresholded! Requested eCprime= ", eCprime, " and lambda1= ", lambda1, "\n")
    perf.c = csnorm:::gfl_perf_iteration(csig, lambda1, lambda2, eCprime)
  }
  if (class(csig)!="CSbdiff") {
    matg = as.data.table(perf.c$mat)
  } else {
    matg = as.data.table(perf.c$mat)
  }
  setkey(matg,bin1,bin2)
  #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()
  return(matg)
}

#' compute BIC for a given value of lambda2, optimizing lambda1 and eCprime (performance iteration, persistent state)
#' @keywords internal
gfl_BIC = function(csig, lambda2, constrained=T, positive=T, fixed=F) {
  #state = perf.c[c("phi.ref","beta")]
  #submat = as.data.table(perf.c$mat)[,.(bin1,bin2,phihat.ref,valuehat=deltahat,ncounts,weight,value=perf.c$delta)]
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
    perf.c = csnorm:::wgfl_signal_BIC(ctsg, dispersion, nperf, nbins, state$GFLState, lambda2,
                                      tol.val, metadata,
                                      state$beta, constrained, fixed)
  } else {
    stopifnot(constrained==T) #for now
    perf.c = csnorm:::wgfl_diff_BIC(ctsg, ctsg.ref, dispersion, nperf, nbins, state$GFLState, lambda2,
                                      tol.val, metadata,
                                      state$phi.ref, state$beta, constrained)
  }
  return(perf.c)
}

#' compute BIC for a given value of lambda1, lambda2 and eCprime (performance iteration, persistent state)
#' @keywords internal
gfl_BIC_fixed = function(csig, lambda1, lambda2, eCprime) {
 stop("signif.threshold=F not implemented!") 
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2 = function(csig, n.SD=1, constrained=T, positive=T, fixed=F, signif.threshold=T, ncores=1) {
  obj = function(x) {
    if (signif.threshold==T) {
      csig@state <<- csnorm:::gfl_BIC(csig, lambda2=10^(x), constrained=constrained, positive=positive, fixed=fixed)
    } else {
      a = csnorm:::gfl_BIC_fixed(csig, 0, lambda2=10^(x), 0)
      if (positive==T) {
        a$eCprime=min(a$beta)
        a$phi=a$beta-min(a$beta)
        a$mat$phi=a$phi
      }
      csig@state <<- a
    }
    #cat("optimize_lambda2: eval at lambda2= ",csig@state$lambda2, " lambda1= ",csig@state$lambda1,
    #    " eCprime= ",csig@state$eCprime," BIC= ",csig@state$BIC, " dof= ",csig@state$dof,"\n")
    return(csig@state$BIC)
  }
  #
  #dt.free = foreach (lam=10^(seq(1.1,1.3,length.out=100)),.combine=rbind) %do% {
  #  obj(log10(lam))
  #  as.data.table(csig@state$mat)[,.(lambda2=lam,lambda1=csig@state$lambda1,eCprime=csig@state$eCprime,bin1,bin2,phihat,phi,beta)]
  #}
  #dt.fix = foreach (lam=10^(seq(log10(12.8),log10(13.3),length.out=20)),.combine=rbind) %dopar% {
  #  csig@state <<- csnorm:::gfl_BIC(csig, lambda2=lam, constrained=constrained, positive=positive, fixed=fixed)
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
  registerDoParallel(cores=ncores)
  for (i in 1:nrounds) {
    minlambda = max(minlambda,lambda2 - range/2)
    maxlambda = min(maxlambda,lambda2 + range/2)
    #cat("round ",i,": min=",minlambda," < lambda2=",lambda2, " < max=",maxlambda," range=",range,"\n")
    dt = rbind(dt,foreach (lam=log10(seq(minlambda, maxlambda, l=npoints)),
                   .combine=rbind) %dopar% data.table(x=lam,y=obj(lam),i=paste(i)))
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
  if (fixed==T && abs(retvals$eCprime)>csig@settings$tol.val) cat("Warning: fixed = T but eCprime != 0\n") #only when signif.threshold==T
  csig@par=modifyList(csig@par,retvals)
  return(csig)
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2_simplified = function(csig, n.SD=1, constrained=T, positive=T, fixed=F, signif.threshold=T, ncores=1) {
  obj = function(x) {
    if (signif.threshold==T) {
      csig@state <<- csnorm:::gfl_BIC(csig, lambda2=10^(x), constrained=constrained, positive=positive, fixed=fixed)
      #cat("optimize at ", 10^(x), " BIC= ", csig@state$BIC,"\n")
    } else {
      a = csnorm:::gfl_BIC_fixed(csig, 0, lambda2=10^(x), 0)
      if (positive==T) {
        a$eCprime=min(a$beta)
        a$phi=a$beta-min(a$beta)
        a$mat$phi=a$phi
      }
      csig@state <<- a
    }
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
  if (fixed==T && abs(retvals$eCprime)>csig@settings$tol.val) cat("Warning: fixed = T but eCprime != 0\n") #only when signif.threshold==T
  csig@par=modifyList(csig@par,retvals)
  return(csig)
}

#' build initial state from phi / delta
#' 
#' @keywords internal
gfl_compute_initial_state = function(csig, diff=F) {
  matg = csig@mat
  if (diff==F) {
    state = list(beta=matg[,phi], lambda1=0.05, eCprime=0, GFLState = list())
  } else {
    state = list(phi.ref=matg[,phi.ref], beta=matg[,delta, GFLState = list()])
  }
  return(state)
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
csnorm_fused_lasso = function(csig, positive, fixed, constrained, verbose=T, signif.threshold=T, ncores=1) {
  n.SD=ifelse(fixed==T,1,0)
  csig = csnorm:::optimize_lambda2_simplified(csig, n.SD=n.SD, constrained=constrained, positive=positive, fixed=fixed,
                                   signif.threshold=signif.threshold, ncores=ncores)
  csig@par$name=csig@cts[,name[1]]
  #matg=csnorm:::gfl_get_matrix(csig, 0, csig@par$lambda2, 0)
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=abs(value-csig@par$eCprime)<=csig@par$lambda1)))
  #matg=csnorm:::gfl_get_matrix(csig, csig@par$lambda1, csig@par$lambda2, csig@par$eCprime)
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  matg = csnorm:::gfl_get_matrix(csig, csig@par$lambda1, csig@par$lambda2, csig@par$eCprime)
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
  mat = csnorm:::get_signal_matrix(cs, resolution, groups=csg@names)
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
  cts=csg@cts[,.(name,bin1,bin2,count,lmu.nosig,weight)]
  csi=new("CSbsig", mat=mat, cts=cts, settings=settings)
  return(csi)
}

#' Build grouped difference matrix using normalization data if available
#' @keywords internal
prepare_difference_estimation = function(cs, csg, resolution, ref, tol.val) {
  csi = csnorm:::prepare_signal_estimation(cs, csg, resolution, tol.val)
  names=csg@names
  mat = foreach(n=names[groupname!=ref,unique(groupname)],.combine=rbind) %do%
    merge(csi@mat[name==n],csi@mat[name==ref,.(bin1,bin2,phi1=phi)],all=T,by=c("bin1","bin2"))
  mat[,c("phi.ref","delta"):=list(phi1,(phi-phi1)/2)]
  mat[,c("phi","phi1"):=NULL]
  cts=csg@cts[name!=ref,.(name,bin1,bin2,count,lmu.nosig,weight)]
  cts.ref=csg@cts[name==ref,.(name,bin1,bin2,count,lmu.nosig,weight)]
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
detect_binless_interactions = function(cs, resolution, group, ncores=1, tol.val=cs@settings$tol.leg, verbose=T, signif.threshold=T){
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
  csi = csnorm:::prepare_signal_estimation(cs, csg, resolution, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  params = foreach(g=groupnames, .combine=rbind) %do% {
    csig = csi
    csig@cts = csi@cts[name==g]
    csig@mat = csi@mat[name==g]
    csig@state = csnorm:::gfl_compute_initial_state(csig, diff=F)
    csnorm:::csnorm_fused_lasso(csig, positive=T, fixed=T, constrained=T, verbose=verbose,
                                signif.threshold=signif.threshold, ncores=ncores)
  }
  #display param info
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],"\n")
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  mat[,c("value","valuehat"):=list(phi,phihat)]
  #report parameters
  csi@par=list(lambda1=params[,lambda1],lambda2=params[,lambda2],name=params[,name])
  #
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_value.png"), width=10, height=8)
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_weight.png"), width=10, height=8)
  #
  if (verbose==T) cat(" Detect patches\n")
  csi@mat = foreach(g=groupnames, .combine=rbind) %do% {
    matg = mat[name==g]
    #cl = csnorm:::build_patch_graph_components(csi@settings$nbins, matg, csi@settings$tol.val)
    #matg[,c("patchno","value"):=list(factor(cl$membership),NULL)]
    matg[,value:=phi]
    matg = csnorm:::detect_binless_patches(matg, csi@settings)
    matg[,value:=NULL]
    matg
  }
  #
  ### store interaction
  #store back
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
detect_binless_differences = function(cs, resolution, group, ref, ncores=1, tol.val=cs@settings$tol.leg, verbose=T, signif.threshold=T){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=", group," and ref=",as.character(ref),"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="CSbdiff", ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  if (is.character(ref)) ref=csg@names[as.character(groupname)==ref,unique(groupname)]
  if (verbose==T) cat("  Prepare for difference estimation\n")
  csi = csnorm:::prepare_difference_estimation(cs, csg, resolution, ref, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  params = foreach(g=groupnames, .combine=rbind) %do% {
    csig = csi
    csig@cts = csi@cts[name==g]
    csig@mat = csi@mat[name==g]
    csig@state = csnorm:::gfl_compute_initial_state(csig, diff=T)
    csnorm:::csnorm_fused_lasso(csig, positive=F, fixed=T, constrained=T, verbose=verbose,
                                signif.threshold=signif.threshold, ncores=ncores)
  }
  #display param info
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],"\n")
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  mat[,c("value","valuehat"):=list(delta,deltahat)]
  #report parameters
  csi@par=list(lambda1=params[,lambda1],lambda2=params[,lambda2],name=params[,name])
  #
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_value.png"), width=10, height=8)
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_weight.png"), width=10, height=8)
  #
  if (verbose==T) cat(" Detect patches\n")
  csi@mat = foreach(g=groupnames, .combine=rbind) %do% {
    matg = mat[name==g]
    #cl = csnorm:::build_patch_graph_components(csi@settings$nbins, matg, csi@settings$tol.val)
    #matg[,c("patchno","value"):=list(factor(cl$membership),NULL)]
    matg[,value:=delta]
    matg = csnorm:::detect_binless_patches(matg, csi@settings)
    matg[,value:=NULL]
    matg
  }
  #store back
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
