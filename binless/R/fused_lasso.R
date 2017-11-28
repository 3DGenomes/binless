#' @include binless.R
NULL

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

#' compute BIC for a given value of lambda2, optimizing lambda1 and eCprime (performance iteration, persistent state)
#' @keywords internal
gfl_BIC = function(csig, lambda2, constrained=T, positive=T, fixed=F, fix.lambda1=F, fix.lambda1.at=NA) {
  stopifnot(fix.lambda1==F || fix.lambda1.at>=0)
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
                                       tol.val, metadata, csig@settings$last.beta,
                                       constrained, fixed, positive, fix.lambda1, fix.lambda1.at)
  } else {
    perf.c = binless:::wgfl_diff_BIC(ctsg, ctsg.ref, dispersion, nperf, nbins, state$GFLState, lambda2,
                                     tol.val, metadata, csig@settings$last.phi.ref, csig@settings$last.beta,
                                     constrained, fix.lambda1, fix.lambda1.at)
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
  # ggplot(dt3[,.SD[1],by=c("lambda2","dset")])+geom_point(aes(lambda2,BIC,colour=converged,group=dset))+facet_grid(dset~., scales="free_y")#+scale_y_log10()
  # ggplot(dt3[,.SD[1],by=c("lambda2","dset")])+geom_line(aes(lambda2,BIC,colour=converged,group=dset))+geom_errorbar(aes(lambda2,ymin=BIC-BIC.sd,ymax=BIC+BIC.sd,colour=converged))+facet_grid(dset~., scales="free_y")#+scale_y_log10()
  # ggplot(dt3[,.SD[BIC==min(BIC)],by=c("dset")])+geom_raster(aes(bin1,bin2,fill=pmin(phihat,3)))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~dset)+coord_fixed()+scale_fill_gradient2()
  # ggplot(dt3[,.SD[abs(lambda2-1)==min(abs(lambda2-1))],by=c("dset")])+geom_raster(aes(bin1,bin2,fill=pmin(phihat,3)))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~dset)+coord_fixed()+scale_fill_gradient2()
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
    l2vals = csig@state$l2vals
    ret = binless:::gfl_BIC(csig, lambda2=10^(x), constrained=constrained, positive=positive, fixed=fixed, fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
    ret$l2vals = rbind(csig@state$l2vals,data.table(lambda2=10^(x),BIC=ret$BIC,BIC.sd=ret$BIC.sd))
    csig@state <<- ret
    #cat("optimize_lambda2: eval at lambda2= ",csig@state$lambda2, " lambda1= ",csig@state$lambda1,
    #    " eCprime= ",csig@state$eCprime," BIC= ",csig@state$BIC, " BIC.sd= ", csig@state$BIC.sd, " dof= ",csig@state$dof,"\n")
    return(csig@state$BIC)
  }
  #first, find rough minimum between 2.5 and 100 by gridding
  minlambda=2.5
  maxlambda=100
  npoints=10
  lvals=10^seq(log10(minlambda),log10(maxlambda),length.out=npoints)
  op = foreach (lam=lvals,.combine=rbind) %do% obj(log10(lam))
  op = copy(csig@state$l2vals)
  setkey(op,lambda2)
  #ggplot(op)+geom_point(aes(lambda2,BIC))+geom_errorbar(aes(lambda2,ymin=BIC-BIC.sd,ymax=BIC+BIC.sd))+scale_x_log10()+scale_y_log10()
  #now find the lowest minimum and extract flanking values
  l2min=op[BIC==min(BIC),min(lambda2)]
  minlambda=op[lambda2<=l2min][.N-1,lambda2]
  if (length(minlambda)==0 || is.na(minlambda)) minlambda=l2min
  maxlambda=op[lambda2>=l2min][2,lambda2]
  if (length(maxlambda)==0 || is.na(maxlambda)) maxlambda=l2min
  #cat("minlambda ",minlambda," maxlambda ",maxlambda,"\n")
  stopifnot(maxlambda>minlambda)
  #optimize there
  op<-optimize(obj, c(log10(minlambda),log10(maxlambda)), tol=0.1)
  lambda2=10^op$minimum
  if (n.SD>0) {
    #finally, find minimum + SD
    optBIC=csig@state$BIC+n.SD*csig@state$BIC.sd
    minlambda=lambda2
    setkey(csig@state$l2vals,lambda2)
    maxlambda=csig@state$l2vals[lambda2>minlambda&BIC>optBIC,lambda2[1]]
    if (length(maxlambda)==0 || is.na(maxlambda)) maxlambda=100
    obj2 = function(x) {
      a=obj(x)
      return(a+2*abs(optBIC-a))
    }
    if (minlambda==maxlambda) {
      lambda2=minlambda
    } else {
      op<-optimize(obj2, c(log10(minlambda),log10(maxlambda)), tol=0.1)
      lambda2=10^op$minimum
    }
  }
  #finish
  if (lambda2==maxlambda) cat("   Warning: lambda2 hit upper boundary.\n")
  if (lambda2==minlambda) cat("   Warning: lambda2 hit lower boundary.\n")
  obj(log10(lambda2))
  retvals = as.list(csig@state)[c("lambda2","lambda1","eCprime","BIC","BIC.sd","dof")]
  if (fixed==T && abs(retvals$eCprime)>csig@settings$tol.val) cat("Warning: fixed = T but eCprime != 0\n") #only when fix.lambda1==F
  csig@par=modifyList(csig@par,retvals)
  csig@state$l2vals=NULL
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
    #first we optimize lambda2 without constraints nor soft-thresholding (setting eCprime=0)
    csig = binless:::optimize_lambda2_smooth(csig, n.SD=n.SD, constrained=F, positive=F, fixed=T,
                                             fix.lambda1=T, fix.lambda1.at=0)
    #now for that optimum, report the best lambda1 (if optimized) and eCprime
    csig = binless:::evaluate_at_lambda2(csig, csig@par$lambda2, constrained=constrained, positive=positive, fixed=fixed,
                                         fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
  } else {
    stopifnot(fix.lambda2.at>0)
    csig = binless:::evaluate_at_lambda2(csig, fix.lambda2.at, constrained=constrained, positive=positive, fixed=fixed,
                                         fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
  }
  #in case we model a constraint, report unconstrained solutions as well
  if (constrained == T) {
    lambda1.constr = csig@par$lambda1
    eCprime.constr = csig@par$eCprime
    BIC.constr = csig@par$BIC
    BIC.sd.constr = csig@par$BIC.sd
    mat.constr = csig@state$mat
    csig = binless:::evaluate_at_lambda2(csig, csig@par$lambda2, constrained=F, positive=positive, fixed=fixed,
                                         fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at)
    csig@par[c("lambda1","lambda1.unconstr","eCprime","eCprime.unconstr","BIC","BIC.unconstr","BIC.sd","BIC.sd.unconstr")]=list(
      lambda1.constr,csig@par$lambda1,eCprime.constr,csig@par$eCprime,BIC.constr,csig@par$BIC,BIC.sd.constr,csig@par$BIC.sd)
    mat.constr$phi.unconstr = csig@state$mat$phi
    csig@state$mat = mat.constr
  }
  csig@par$name=csig@cts[,name[1]]
  matg = as.data.table(csig@state$mat)
  setkey(matg,bin1,bin2)
  matg[,name:=csig@par$name]
  params = as.data.table(csig@par)
  params[,c("state","mat"):=list(list(csig@state),list(matg))]
  return(params)
}

