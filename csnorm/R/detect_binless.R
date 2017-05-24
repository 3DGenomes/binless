#' @include csnorm.R
NULL

#' give a patch ID to each patch in a binless matrix, and report local extrema
#' @keywords internal
build_patch_graph = function(mat, trails, tol.value=1e-3) {
  g=trails$graph
  #remove edges that connect different values
  V(g)$value = mat[,value]
  values_by_edge = t(matrix(get.vertex.attribute(g,"value",index=get.edges(g,E(g))), ncol=2))
  to_rm = E(g)[abs(values_by_edge[1,]-values_by_edge[2,])>tol.value]
  g2 = delete_edges(g, to_rm)
  #plot(g2, vertex.color=factor(V(g2)$value), layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #
  #get connected components
  cl=clusters(g2)
  #V(g)$patch=cl$membership
  #plot(g, vertex.color=factor(V(g)$value), vertex.label=V(g)$patch, layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #V(g2)$patch=cl$membership
  #plot(g2, vertex.color=factor(V(g2)$value), vertex.label=V(g2)$patch, layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #patches=foreach (id=1:cl$no) %do% V(g)[cl$membership==id]
  #plot(g, vertex.color=factor(V(g)$value), vertex.label=V(g)$patch, edge.arrow.size=0.3,
  #     mark.groups=patches, layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  return(list(g.patches=g2,components=cl))
}

#' give a patch ID to each patch in a binless matrix, and report local extrema
#' @keywords internal
detect_binless_patches = function(mat, trails, tol.value=1e-3) {
  if (mat[,uniqueN(name)]>1)
    return(foreach (n=mat[,unique(name)],.combine=rbind) %do%
             csnorm:::detect_binless_patches(mat[name==n], trails, tol.value=tol.value))
  #
  stuff = build_patch_graph(mat, trails, tol.value=tol.value)
  g=trails$graph
  V(g)$value = mat[,value]
  cl=stuff$components
  #create fused graph
  g3 = contract(g, cl$membership, vertex.attr.comb = list(value="min", weight="sum")) %>% simplify()
  #plot(g3, vertex.color=factor(V(g3)$value), vertex.label=V(g3)$weight, vertex.size=10*log(1+V(g3)$weight), layout=layout_as_tree)
  #
  #Add directions from high to low values
  g4 = as.directed(g3)
  values_by_edge = t(matrix(get.vertex.attribute(g4,"value",index=get.edges(g4,E(g4))), ncol=2))
  to_rm = E(g4)[values_by_edge[1,]<values_by_edge[2,]]
  g4 = delete_edges(g4, to_rm)
  #plot(g4, vertex.color=factor(V(g4)$value), vertex.label=V(g4)$weight,
  #     vertex.size=5*V(g4)$weight, edge.arrow.size=0.5, layout=layout_nicely)
  #
  #report in and outdegree to original graph
  g5 = as.directed(g)
  V(g5)$indegree=degree(g4, v=cl$membership, mode="in")
  V(g5)$outdegree=degree(g4, v=cl$membership, mode="out")
  values_by_edge = t(matrix(get.vertex.attribute(g5,"value",index=get.edges(g5,E(g5))), ncol=2))
  to_rm = E(g5)[values_by_edge[1,]<=values_by_edge[2,]]
  g5 = delete_edges(g5, to_rm)
  #plot(g5, vertex.color=V(g5)$indegree==0, vertex.label=NA, edge.arrow.size=0.3, vertex.size=3,
  #     mark.groups=patches, layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #plot(g5, vertex.color=V(g5)$outdegree==0, vertex.label=NA, edge.arrow.size=0.3,
  #     mark.groups=patches, layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #source_or_sink=ifelse(V(g5)$indegree==0,1,ifelse(V(g5)$outdegree==0,2,0))
  #plot(g5, vertex.color=source_or_sink, vertex.label=NA, edge.arrow.size=0.01, vertex.size=3,
  #     layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  stopifnot(length(V(g5))==mat[,.N])
  mat[,is.minimum:=V(g5)$outdegree==0]
  mat[,is.maximum:=V(g5)$indegree==0]
  mat[,patchno:=factor(cl$membership)]
  return(mat)
}

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

#' compute trail information for the gfl package, for a triangle trid with nrow rows
#' @keywords internal
gfl_compute_trails = function(nrow) {
  chain = csnorm:::boost_triangle_grid_chain(nrow)
  trails = csnorm:::boost_chains_to_trails(chain)
  stopifnot(uniqueN(trails$trails)==nrow*(nrow+1)/2)
  #store bin graph
  trails$graph = csnorm:::compute_2d_connectivity_graph(nrow)
  #plot(trails$graph, vertex.color=factor(V(g)$value), vertex.label=V(g)$count,
  #     layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  trails$nrow=nrow
  return(trails)
}

#' dispatch for fused lasso perf iteration: warm/cold signal/difference
#' @keywords internal
gfl_perf_iteration = function(csig, lambda2) {
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
    if (length(state)==0) {
      perf.c = csnorm:::wgfl_signal_perf(ctsg, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                         trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.val/20, diag.rm)
    } else {
      perf.c = csnorm:::wgfl_signal_perf_warm(ctsg, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                              trails$breakpoints, lambda2, state$alpha, inflate, maxsteps, tol.val/20, diag.rm,
                                              state$z, state$u, state$phi)
    }
  } else {
    if (length(state)==0) {
      perf.c = csnorm:::wgfl_diff_perf(ctsg, ctsg.ref, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                       trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.val/20, diag.rm)
    } else {
      perf.c = csnorm:::wgfl_diff_perf_warm(ctsg, ctsg.ref, dispersion, nperf, nbins, trails$ntrails, trails$trails,
                                            trails$breakpoints, lambda2, state$alpha, inflate, maxsteps, tol.val/20, diag.rm,
                                            state$z, state$u, state$phi.ref, state$delta)
    }
  }
  return(perf.c)
}

#' compute sparse fused lasso matrix for a given value of lambda1, lambda2 and eCprime (performance iteration)
#' @keywords internal
gfl_get_matrix = function(csig, lambda1, lambda2, eCprime) {
  #assume lambda1=0 and compute the fused lasso solution, centered around eCprime
  perf.c = csnorm:::gfl_perf_iteration(csig, lambda2)
  #cat("GFL: initial alpha=",alpha,"final alpha=",perf.c$alpha,"nsteps=",perf.c$nsteps,"\n")
  if (class(csig)!="CSbdiff") {
    matg = as.data.table(perf.c$mat)[,.(bin1,bin2,valuehat=phihat,ncounts,weight,value=perf.c$phi-eCprime,diag.idx)]
  } else {
    matg = as.data.table(perf.c$mat)[,.(bin1,bin2,phihat.ref,valuehat=deltahat,ncounts,weight,value=perf.c$delta-eCprime,diag.idx)]
  }
  #now soft-threshold the shifted value around eCprime
  matg[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
  #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()
  return(matg)
}

#' compute BIC for a given value of lambda1, lambda2 and eCprime (performance iteration, persistent state)
#' @keywords internal
gfl_BIC = function(csig, lambda1, lambda2, eCprime) {
  #get value with lambda1 set to zero to avoid round-off errors in degrees of freedom
  perf.c = csnorm:::gfl_perf_iteration(csig, lambda2)
  if (class(csig)!="CSbdiff") {
    state = perf.c[c("z","u","phi","alpha")]
    submat = as.data.table(perf.c$mat)[,.(bin1,bin2,valuehat=phihat,ncounts,weight,value=perf.c$phi-eCprime)]
  } else {
    state = perf.c[c("z","u","phi.ref","delta","alpha")]
    submat = as.data.table(perf.c$mat)[,.(bin1,bin2,phihat.ref,valuehat=deltahat,ncounts,weight,value=perf.c$delta-eCprime)]
  }
  #get the number of patches and deduce degrees of freedom
  tol.value=csig@settings$tol.val
  cl = csnorm:::build_patch_graph(submat, csig@trails, tol.value=tol.value)$components
  submat[,patchno:=cl$membership]
  dof = submat[abs(value)+tol.value>lambda1,uniqueN(patchno)] #sparse fused lasso
  state$dof=dof
  stopifnot(dof<=cl$no)
  #now soft-threshold the value around eCprime
  submat[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
  #compute BIC
  BIC = submat[,sum(weight*((valuehat-(value+eCprime))^2))+log(sum(ncounts))*dof]
  #compute mallow's Cp
  #Cp = submat[,sum(weight*((valuehat-(value+eCprime))^2 - 1))]+2*dof
  state$BIC=BIC
  state$mat=submat
  return(state)
}

#' cross-validate lambda1 and eCprime, assuming positive=T and holding eCprime=lambda1+min(value)
#' 
#' @keywords internal
optimize_lambda1_eCprime = function(matg, csig, lambda1.min=0.05, constrained=T) {
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, csig@trails, tol.value=csig@settings$tol.val)$components
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1],size=.N),by=patchno][,.(N=.N,size=size[1]),keyby=value]
  if (patches[,.N]==1) {#matrix is just one big bin
    csig@par=modifyList(csig@par,list(eCprime=patches[,value],lambda1=lambda1.min,
                                      BIC=matg[,sum(weight*((valuehat-(value+patches[,value]))^2))],dof=0))
    return(csig)
  }
  stopifnot(patches[,sum(N)]==cl$no)
  minval=patches[,min(value)]
  maxval=patches[,max(value)]
  valrange=maxval-minval
  #
  if (constrained==T) { #some patches must be zeroed to avoid degeneracy with diagonal decay fit
    forbidden.vals = matg[,min(value.ori),by=diag.idx][,unique(V1)]
    lambda1.min = max(lambda1.min, (max(forbidden.vals)-minval)/2 + csig@settings$tol.val)
  } else {
    forbidden.vals = c()
  }
  #
  obj = function(lambda1) {
    if (valrange<2*lambda1+csig@settings$tol.val) {
      eCprime=(maxval+minval)/2
    } else {
      eCprime=lambda1+minval-csig@settings$tol.val
    }
    if (constrained==T & any(abs(eCprime-forbidden.vals)>lambda1+csig@settings$tol.val))
      return(data.table(eCprime=0,lambda1=lambda1,BIC=.Machine$double.xmax,dof=NA))
    matg[,value:=value.ori-eCprime]
    dof = matg[abs(value)>lambda1,uniqueN(patchno)] #sparse fused lasso
    stopifnot(dof<=cl$no)
    #now soft-threshold the value around eCprime
    matg[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
    #compute BIC
    #BIC = matg[,sum(weight*((valuehat-(value+eCprime))^2))+log(sum(ncounts))*dof-2*.N*log(lambda1)]
    BIC = matg[,sum(weight*((valuehat-(value+eCprime))^2))+log(sum(ncounts))*dof]#-9*log(lambda1)+5*lambda1]
    data.table(eCprime=eCprime,lambda1=lambda1,BIC=BIC,dof=dof)
  }
  #dt=foreach(lambda1=seq(lambda1.min,valrange,length.out=50), .combine=rbind) %do% obj(lambda1)
  #ggplot(dt)+geom_point(aes(lambda1,BIC))+geom_line(aes(lambda1,BIC))
  if (valrange <= 2*lambda1.min) {
    op=list(minimum=lambda1.min)
  } else {
    op=optimize(function(x){obj(10^(x))[,BIC]}, c(log10(max(lambda1.min,csig@settings$tol.val/2)),log10(valrange)),
              tol=csig@settings$tol.val)
  }
  csig@par=modifyList(csig@par,as.list(obj(10^(op$minimum))))
  #patches[,removed:=abs(value-values$eCprime)<=values$lambda1]
  return(csig)
}

#' cross-validate lambda1 and assume eCprime=0
#' 
#' @keywords internal
optimize_lambda1_only = function(matg, csig, lambda1.min=0.05, positive=F, constrained=T) {
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, csig@trails, tol.value=csig@settings$tol.val)$components
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
    csig@state <<- csnorm:::gfl_BIC(csig, lambda1=0, lambda2=10^(x), eCprime=0)
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
  csig@par$lambda2=lambda2
  return(csig)
}

#' build initial state from phi / delta
#' 
#' @keywords internal
gfl_compute_initial_state = function(csig, diff=F, init.alpha=5) {
  matg = csig@mat
  trails = csig@trails
  if (diff==F) {
    state = list(phi=matg[,phi], u=rep(0,length(trails$trails)), z=matg[trails$trails+1,phi], alpha=init.alpha)
  } else {
    state = list(phi.ref=matg[,phi.ref], delta=matg[,delta], u=rep(0,length(trails$trails)),
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
csnorm_fused_lasso = function(csig, positive, fixed, constrained, verbose=T, ctsg.ref=NULL) {
  csig = csnorm:::optimize_lambda2(csig)
  #compute values for lambda1=0 and eCprime=0
  matg = csnorm:::gfl_get_matrix(csig, 0, csig@par$lambda2, 0)
  matg[,value.ori:=value]
  #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()
  #get best lambda1 and set eCprime to lower bound
  if (fixed==F) {
    csig = csnorm:::optimize_lambda1_eCprime(matg, csig, constrained=constrained)
  } else {
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
prepare_signal_matrix = function(cs, names, resolution) {
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
  trails = csnorm:::gfl_compute_trails(mat[,nlevels(bin1)])
  stopifnot(all(mat[,.N,by=name]$N==mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  stopifnot(all(length(V(trails$graph))==mat[,.N,by=name]$N))
  csi=new("CSbsig", mat=mat, trails=trails, cts=data.table())
  return(csi)
}

#' Build grouped difference matrix using normalization data if available
#' @keywords internal
prepare_difference_matrix = function(cs, names, resolution, ref) {
  csi = csnorm:::prepare_signal_matrix(cs, names, resolution)
  mat = foreach(n=names[groupname!=ref,unique(groupname)],.combine=rbind) %do%
    merge(csi@mat[name==n],csi@mat[name==ref,.(bin1,bin2,phi1=phi)],all=T,by=c("bin1","bin2"))
  mat[,c("phi.ref","delta"):=list(phi1,(phi-phi1)/2)]
  mat[,c("phi","phi1"):=NULL]
  csi=new("CSbdiff", mat=mat, trails=csi@trails, cts=data.table(), cts.ref=data.table(), ref=as.character(ref))
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
  csi = csnorm:::prepare_signal_matrix(cs, csg@names, resolution)
  csi@settings$diag.rm = ceiling(csg@par$dmin/resolution)
  csi@settings$nbins = csg@par$nbins
  csi@settings$dispersion = csg@par$alpha
  csi@settings$tol.val = tol.val
  csi@settings$inflate=2
  csi@settings$nperf=1000
  csi@settings$maxsteps=100000
  csi@cts=csg@cts[,.(name,bin1,bin2,count,lmu.nosig,weight)]
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
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
          " eCprime=",params[i,eCprime],"\n")
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
  csi@mat = csnorm:::detect_binless_patches(mat, csi@trails, tol.value=csi@settings$tol.val)
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
  csi = csnorm:::prepare_difference_matrix(cs, csg@names, resolution, ref)
  csi@settings$diag.rm = ceiling(csg@par$dmin/resolution)
  csi@settings$nbins = csg@par$nbins
  csi@settings$dispersion = csg@par$alpha
  csi@settings$tol.val = tol.val
  csi@settings$inflate=2
  csi@settings$nperf=1000
  csi@settings$maxsteps=100000
  csi@cts=csg@cts[name!=ref,.(name,bin1,bin2,count,lmu.nosig,weight)]
  csi@cts.ref=csg@cts[name==ref,.(name,bin1,bin2,count,lmu.nosig,weight)]
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
      cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
          " eCprime=",params[i,eCprime],"\n")
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
  csi@mat = csnorm:::detect_binless_patches(mat, csi@trails, tol.value=csi@settings$tol.val)
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
