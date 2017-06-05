#' @include csnorm.R
NULL

#' connectivity on a triangle
#' @keywords internal
compute_2d_connectivity_graph = function(nbins, start=1) {
  stopifnot(start>=1,nbins>=2)
  if (nbins==2) {
    ret=graph(c(start,start+1,start+1,start+2), directed=F)
  } else {
    upper.row = c(c(start,start+1),
                  c(sapply((start+1):(start+nbins-2),function(x){c(x,x+1,x,x+nbins-1)},simplify=T)),
                  c(start+nbins-1,start+2*(nbins-1)))
    lower.tri = compute_2d_connectivity_graph(nbins-1,start=start+nbins)
    ret = add.edges(lower.tri,upper.row)
  }
  #plot(ret)
  return(ret)
}

#' compute trail information for the gfl package, for a triangle trid with nbins rows
#' @keywords internal
gfl_compute_trails = function(nbins) {
  chain = csnorm:::boost_triangle_grid_chain(nbins)
  trails = csnorm:::boost_chains_to_trails(chain)
  stopifnot(uniqueN(trails$trails)==nbins*(nbins+1)/2)
  #store bin graph
  trails$graph = csnorm:::compute_2d_connectivity_graph(nbins)
  #plot(trails$graph, vertex.color=factor(V(g)$value), vertex.label=V(g)$count,
  #     layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  return(trails)
}

#' dispatch for fused lasso perf iteration: warm/cold signal/difference
#' @keywords internal
gfl_perf_iteration = function(csig, lambda1, lambda2, eCprime) {
  ctsg=csig@cts
  nbins=csig@settings$nbins
  dispersion=csig@settings$dispersion
  diag.rm=csig@settings$diag.rm
  trails=csig@trails
  tol.val=csig@settings$tol.val
  state=csig@state
  stopifnot(length(state)>0) #warm start only
  if (class(csig)=="CSbdiff") ctsg.ref=csig@cts.ref else ctsg.ref=NULL
  inflate=csig@settings$inflate
  nperf=csig@settings$nperf
  maxsteps=csig@settings$maxsteps
  if (is.null(ctsg.ref)) {
    perf.c = csnorm:::wgfl_signal_perf_warm(ctsg, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                              trails$breakpoints, lambda1, lambda2, eCprime,
                                              state$alpha, inflate, maxsteps, tol.val/20, diag.rm,
                                              state$z, state$u, state$beta)
  } else {
    stopifnot(lambda1==0 & eCprime==0)
    perf.c = csnorm:::wgfl_diff_perf_warm(ctsg, ctsg.ref, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                            trails$breakpoints, lambda2, state$alpha, inflate, maxsteps, tol.val/20, diag.rm,
                                            state$z, state$u, state$phi.ref, state$delta)
  }
  return(perf.c)
}

#' compute sparse fused lasso matrix for a given value of lambda1, lambda2 and eCprime (performance iteration)
#' @keywords internal
gfl_get_matrix = function(csig, lambda1, lambda2, eCprime) {
  #assume lambda1=0 and compute the fused lasso solution, centered around eCprime
  perf.c = csnorm:::gfl_perf_iteration(csig, lambda1, lambda2, eCprime)
  #cat("GFL: initial alpha=",alpha,"final alpha=",perf.c$alpha,"nsteps=",perf.c$nsteps,"\n")
  if (class(csig)!="CSbdiff") {
    matg = as.data.table(perf.c$mat)[,.(bin1,bin2,valuehat=phihat,ncounts,weight,value=perf.c$phi,diag.idx)]
  } else {
    matg = as.data.table(perf.c$mat)[,.(bin1,bin2,phihat.ref,valuehat=deltahat,ncounts,weight,value=perf.c$delta,diag.idx)]
  }
  setkey(matg,bin1,bin2)
  #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()
  return(matg)
}

#' compute BIC for a given value of lambda1, lambda2 and eCprime (performance iteration, persistent state)
#' @keywords internal
gfl_BIC = function(csig, lambda2, lambda1.min=0.05, percent.closest=50) {
  stopifnot(class(csig)!="CSbdiff")
  #state = perf.c[c("z","u","phi.ref","beta","alpha")]
  #submat = as.data.table(perf.c$mat)[,.(bin1,bin2,phihat.ref,valuehat=deltahat,ncounts,weight,value=perf.c$delta)]
  #
  ctsg=csig@cts
  nbins=csig@settings$nbins
  dispersion=csig@settings$dispersion
  diag.rm=csig@settings$diag.rm
  trails=csig@trails
  tol.val=csig@settings$tol.val
  state=csig@state
  stopifnot(length(state)>0) #warm start only
  if (class(csig)=="CSbdiff") ctsg.ref=csig@cts.ref else ctsg.ref=NULL
  inflate=csig@settings$inflate
  nperf=csig@settings$nperf
  maxsteps=csig@settings$maxsteps
  perf.c = csnorm:::wgfl_signal_BIC(ctsg, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                    trails$breakpoints, lambda2,
                                    state$alpha, inflate, maxsteps, tol.val, diag.rm,
                                    state$z, state$u, state$beta, lambda1.min, percent.closest)
  return(perf.c)
}

#' cross-validate lambda1 and assume eCprime=0
#' 
#' @keywords internal
optimize_lambda1_only = function(matg, csig, lambda1.min=0.05, positive=F, constrained=T) {
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::boost_build_patch_graph_components(csig@settings$nbins, matg, csig@settings$tol.val)
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1],size=.N),by=patchno][,.(N=.N,size=size[1]),keyby=value]
  stopifnot(patches[,sum(N)]==cl$no)
  if (positive==T) {
    minval = max(patches[,abs(min(value))], lambda1.min)
  } else {
    minval = max(patches[,min(abs(value))], lambda1.min)
  }
  maxval=patches[,max(abs(value))]
  #
  if (constrained==T) { #some patches must be zeroed to avoid degeneracy with diagonal decay fit
    if (positive==T) {
      forbidden.vals = matg[,min(value.ori),by=diag.idx][,unique(V1)]
    } else {
      forbidden.vals = matg[,.(diff=max(value.ori)-min(value.ori),val=min(value.ori)+tol.val/2),by=diag.idx][
        diff<=csig@settings$tol.val,unique(val)]
    }
    minval = max(minval, abs(forbidden.vals))
  } else {
    forbidden.vals = c()
  }
  #
  if (patches[,.N]==1) if (abs(patches[,value])<=minval) {#matrix is just one big bin
    csig@par=modifyList(csig@par,list(eCprime=0, lambda1=minval,
                                      BIC=matg[,sum(weight*(valuehat-value)^2)],dof=0))
    return(csig)
  }
  #
  obj = function(lambda1) {
    dof = matg[abs(value.ori)>lambda1,uniqueN(patchno)] #sparse fused lasso
    stopifnot(dof<=cl$no)
    #now soft-threshold the value around zero
    matg[,value:=sign(value.ori)*pmax(abs(value.ori)-lambda1, 0)]
    #compute BIC
    #BIC = matg[,sum(weight*((valuehat-(value+eCprime))^2))+log(sum(ncounts))*dof-2*.N*log(lambda1)]
    BIC = matg[,sum(weight*((valuehat-value)^2))+log(sum(ncounts))*dof]#-9*log(lambda1)+5*lambda1]
    data.table(eCprime=0,lambda1=lambda1,BIC=BIC,dof=dof)
  }
  #dt=foreach(lambda1=seq(lambda1.min,maxval,length.out=50), .combine=rbind) %do% obj(lambda1)
  #ggplot(dt)+geom_point(aes(lambda1,BIC))+geom_line(aes(lambda1,BIC))
  op=optimize(function(x){obj(10^(x))[,BIC]}, c(log10(max(minval,csig@settings$tol.val/2)),log10(maxval)),
              tol=csig@settings$tol.val)
  csig@par=modifyList(csig@par,as.list(obj(10^(op$minimum))))
  return(csig)
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2 = function(csig, minlambda=0.1, maxlambda=100) {
  obj = function(x) {
    csig@state <<- csnorm:::gfl_BIC(csig, lambda2=10^(x))
    return(csig@state$BIC)
  }
  #ctsg=copy(ctsg.old)
  #dt.old = foreach (lam=seq(log10(minlambda),log10(maxlambda),l=100),.combine=rbind) %dopar% data.table(x=lam,y=obj(lam))
  #ctsg=copy(ctsg.new)
  #dt.new = foreach (lam=seq(log10(minlambda),log10(maxlambda),l=100),.combine=rbind) %dopar% data.table(x=lam,y=obj(lam))
  #ggplot(rbind(dt.old[,.(x,y,ori="old")],dt.new[,.(x,y,ori="new")]))+geom_line(aes(10^(x),y,colour=ori))#+scale_y_log10()
  op=optimize(obj, c(log10(minlambda),log10(maxlambda)), tol=csig@settings$tol.val)
  lambda2=10^op$minimum
  if (lambda2==maxlambda) cat("   Warning: lambda2 hit upper boundary.\n")
  if (lambda2 <= csig@settings$tol.val*10) {
    cat("   Warning: lambda2 too close to lower boundary.")
    lambda2=0
  }
  retvals = as.list(csnorm:::gfl_BIC(csig, lambda2))[c("lambda2","lambda1","eCprime","BIC","dof")]
  csig@par=modifyList(csig@par,retvals)
  return(csig)
}

#' build initial state from phi / delta
#' 
#' @keywords internal
gfl_compute_initial_state = function(csig, diff=F, init.alpha=5) {
  matg = csig@mat
  trails = csig@trails
  if (diff==F) {
    state = list(beta=matg[,phi], u=rep(0,length(trails$trails)), z=matg[trails$trails+1,phi], alpha=init.alpha)
  } else {
    state = list(phi.ref=matg[,phi.ref], beta=matg[,delta], u=rep(0,length(trails$trails)),
                 z=matg[trails$trails+1,delta], alpha=init.alpha)
  }
  return(state)
}

#' run fused lasso on one dataset contained in matg, fusing 'valuehat' into
#' 'value'
#' 
#' @param matg a data.table containing one dataset
#' @param trails the trails list at that resolution
#' @param positive boolean. Constrain eCprime in order to force 'value' to be
#'   positive?
#' @param fixed boolean. Set eCprime=0 throughout ?
#' @param constrained boolean. Constrain lambda1 so that any diagonal that contains
#'   only one big patch be forced to have 'value'=0 ?
#' @param tol.val numeric. Convergence tolerance on fused value
#' @param verbose boolean (default TRUE).
#' @param init.state initial state for lasso, for speedup purposes. Set to NULL if unknown.
#' @param ctsg.ref if provided, compute fused lasso on differences wrt this dataset
#'   
#'   finds optimal lambda1, lambda2 and eC using BIC.
#' @keywords internal
csnorm_fused_lasso = function(csig, positive, fixed, constrained, verbose=T, ctsg.ref=NULL, lambda1.min=0.05) {
  stopifnot(positive==T & fixed==F)
  csig = csnorm:::optimize_lambda2(csig)
  #compute values for lambda1=0 and eCprime=0
  matg = csnorm:::gfl_get_matrix(csig, 0, csig@par$lambda2, 0)
  matg[,value.ori:=value]
  #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()
  #get best lambda1 and set eCprime to lower bound
  if (fixed==T) {
    csig = csnorm:::optimize_lambda1_only(matg, csig, constrained=constrained, positive=positive)
  }
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

#' Build grouped signal matrix using normalization data if available
#' @keywords internal
prepare_signal_matrix = function(cs, csg, resolution, tol.val) {
  names=csg@names
  #build signal matrix
  if (cs@par$signal[,.N]==0) {
    sbins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
    signal.bins=unique(cut(c(sbins,head(sbins,n=-1)+resolution/2), sbins,
                           ordered_result=T, right=F, include.lowest=T,dig.lab=12))
    mat=CJ(name=names[,groupname],bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
    mat[,phi:=0]
  } else {
    #report phi values for each group
    if (resolution!=cs@settings$base.res) {
      refmat=names[cs@par$signal[,.(name,bin1,bin2,phi)]][
        !is.na(groupname),.(name=groupname,refbin1=bin1,refbin2=bin2,phi)][
          ,.(phi=mean(phi)),by=c("name","refbin1","refbin2")]
      #merge signal to new binning
      sbins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
      pos=head(sbins,n=-1)+resolution/2
      bins=unique(data.table(refbin=cut(pos, cs@settings$sbins,
                                        ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                             bin=cut(pos, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)))
      stopifnot(bins[,.N]==length(sbins)-1)
      mat=CJ(name=names[,groupname],bin1=bins[,unique(bin)],bin2=bins[,unique(bin)],
             sorted=F,unique=F)[bin2>=bin1]
      mat=merge(mat,bins,by.x="bin1",by.y="bin",all.y=T)
      mat=merge(mat,bins,by.x="bin2",by.y="bin",all.y=T, suffixes=c("1","2"))
      mat=merge(mat,refmat,by=c("name","refbin1","refbin2"),all.x=T)
      mat[is.na(phi),phi:=0]
      mat=mat[,.(name,bin1,bin2,phi)]
    } else {
      mat=names[cs@par$signal[,.(name,bin1,bin2,phi)]][!is.na(groupname),.(name=groupname,bin1,bin2,phi)]
    }
    mat=mat[,.(phi=mean(phi)),keyby=c("name","bin1","bin2")]
  }
  #ggplot(cs@par$signal)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
  #
  #add trail information
  trails = csnorm:::gfl_compute_trails(csg@par$nbins)
  stopifnot(all(mat[,.N,by=name]$N==mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  stopifnot(all(length(V(trails$graph))==mat[,.N,by=name]$N))
  #add other settings
  settings=list(diag.rm = ceiling(csg@par$dmin/resolution),
                nbins = csg@par$nbins,
                dispersion = csg@par$alpha,
                tol.val = tol.val,
                inflate=2,
                nperf=100,
                maxsteps=100000)
  cts=csg@cts[,.(name,bin1,bin2,count,lmu.nosig,weight)]
  csi=new("CSbsig", mat=mat, trails=trails, cts=cts, settings=settings)
  return(csi)
}

#' Build grouped difference matrix using normalization data if available
#' @keywords internal
prepare_difference_matrix = function(cs, csg, resolution, ref, tol.val) {
  csi = csnorm:::prepare_signal_matrix(cs, csg, resolution, tol.val)
  names=csg@names
  mat = foreach(n=names[groupname!=ref,unique(groupname)],.combine=rbind) %do%
    merge(csi@mat[name==n],csi@mat[name==ref,.(bin1,bin2,phi1=phi)],all=T,by=c("bin1","bin2"))
  mat[,c("phi.ref","delta"):=list(phi1,(phi-phi1)/2)]
  mat[,c("phi","phi1"):=NULL]
  cts=csg@cts[name!=ref,.(name,bin1,bin2,count,lmu.nosig,weight)]
  cts.ref=csg@cts[name==ref,.(name,bin1,bin2,count,lmu.nosig,weight)]
  csi=new("CSbdiff", mat=mat, trails=csi@trails, cts=cts, cts.ref=cts.ref,
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
detect_binless_interactions = function(cs, resolution, group, ncores=1, tol.val=1e-5, verbose=T){
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
  csi = csnorm:::prepare_signal_matrix(cs, csg, resolution, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  registerDoParallel(cores=ncores)
  params = foreach(g=groupnames, .combine=rbind) %dopar% {
    csig = csi
    csig@cts = csi@cts[name==g]
    csig@mat = csi@mat[name==g]
    csig@state = csnorm:::gfl_compute_initial_state(csig, diff=F, init.alpha=5)
    csnorm:::csnorm_fused_lasso(csig, positive=T, fixed=T, constrained=T, verbose=verbose)
  }
  #display param info
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],"\n")
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  mat[,phi:=value]
  #
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_value.png"), width=10, height=8)
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_weight.png"), width=10, height=8)
  #
  if (verbose==T) cat(" Detect patches\n")
  csi@mat = foreach(g=groupnames, .combine=rbind) %dopar% {
    matg = csi@mat[name==g,.(name,bin1,bin2,value=phi)]
    cl = csnorm:::boost_build_patch_graph_components(csi@settings$nbins, matg, csi@settings$tol.val)
    matg[,patchno:=factor(cl$membership)]
    setnames(matg,"value","phi")
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
detect_binless_differences = function(cs, resolution, group, ref, ncores=1, tol.val=1e-3, verbose=T){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=", group," and ref=",as.character(ref),"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="CSbdiff", ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  if (is.character(ref)) ref=csg@names[as.character(groupname)==ref,unique(groupname)]
  if (verbose==T) cat("  Prepare for difference estimation\n")
  csi = csnorm:::prepare_difference_matrix(cs, csg, resolution, ref, tol.val)
  #
  #perform fused lasso on signal
  if (verbose==T) cat("  Fused lasso\n")
  groupnames=csi@cts[,unique(name)]
  registerDoParallel(cores=ncores)
  params = foreach(g=groupnames, .combine=rbind) %dopar% {
    csig = csi
    csig@cts = csi@cts[name==g]
    csig@mat = csi@mat[name==g]
    csig@state = csnorm:::gfl_compute_initial_state(csig, diff=T, init.alpha=5)
    csnorm:::csnorm_fused_lasso(csig, positive=F, fixed=T, constrained=T, verbose=verbose, ctsg.ref=csig@cts.ref)
  }
  #display param info
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],"\n")
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  mat[,delta:=value]
  #
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_value.png"), width=10, height=8)
  #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
  #ggsave(p,filename = paste0("sig_step_",step,"_weight.png"), width=10, height=8)
  #
  if (verbose==T) cat(" Detect patches\n")
  csi@mat = foreach(g=groupnames, .combine=rbind) %dopar% {
    matg = csi@mat[name==g,.(name,bin1,bin2,value=delta)]
    cl = csnorm:::boost_build_patch_graph_components(csi@settings$nbins, matg, csi@settings$tol.val)
    matg[,patchno:=factor(cl$membership)]
    setnames(matg,"value","delta")
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
plot_binless_matrix = function(mat, minima=F, scale=T) {
  resolution=mat[bin1==bin1[1]&name==name[1],begin2[2]-begin2[1]]
  a=mat[is.maximum==T]
  a=a[,.SD[,.(begin1=c(begin1,begin1,end1,end1)-resolution/2, begin2=c(begin2,end2,begin2,end2)-resolution/2,
              patchno, value)][chull(begin1,begin2)], by=c("patchno","name")]
  if (minima==T) {
    b=mat[is.minimum==T]
    b=b[,.SD[,.(begin1=c(begin1,begin1,end1,end1)-resolution/2, begin2=c(begin2,end2,begin2,end2)-resolution/2,
                patchno, value)][chull(begin1,begin2)], by=c("patchno","name")]
    p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=-(value)))+
      geom_raster(aes(begin2,begin1,fill=-(value)))+
      scale_fill_gradient2()+facet_wrap(~name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="blue",fill=NA,data=b)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="red",fill=NA,data=a)
  } else {
    p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=-value))+
      geom_raster(aes(begin2,begin1,fill=-value))+
      scale_fill_gradient2()+facet_wrap(~name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="black",fill=NA,data=a)
  }
  if (scale!=T) p=p+guides(fill=F)
  print(p)
}
