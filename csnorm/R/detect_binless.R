#' @include csnorm.R
NULL

#' internal use
#' @keywords internal
gfl_triangle_grid_chain = function(nrow) {
  #rows of consecutive numbers
  ntotal = nrow*(nrow+1)/2-1
  l=nrow
  chains=list()
  current=c(0)
  for (i in 1:ntotal) {
    if (length(current)==l) {
      chains=c(chains,list(current))
      current=c(i)
      l=l-1
    } else {
      current=c(current,i)
    }
  }
  #columns with Ui+1 = Ui + (N-i) with U1 from 2 to nrow
  for (U1 in 2:nrow) {
    Ui=U1
    current=c(Ui-1)
    for (i in 1:(U1-1)) {
      Uip1=Ui+nrow-i
      current=c(current,Uip1-1)
      Ui=Uip1
    }
    chains=c(chains,list(current))
  }
  return(chains)
}

#' internal use
#' @keywords internal
gfl_chains_to_trails = function(chains) {
  trails=c()
  breakpoints=c()
  for (nodes in chains) {
    if (length(trails)>0) breakpoints=c(breakpoints,length(trails))
    trails=c(trails,nodes)
  }
  if (length(trails)>0) breakpoints=c(breakpoints,length(trails))
  return(list(ntrails=length(breakpoints),trails=trails,breakpoints=breakpoints))
}

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
  chain = csnorm:::gfl_triangle_grid_chain(nrow)
  trails = csnorm:::gfl_chains_to_trails(chain)
  stopifnot(uniqueN(trails$trails)==nrow*(nrow+1)/2)
  #store bin graph
  trails$graph = csnorm:::compute_2d_connectivity_graph(nrow)
  #plot(trails$graph, vertex.color=factor(V(g)$value), vertex.label=V(g)$count,
  #     layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  return(trails)
}

#' compute sparse fused lasso coefficients for a given value of lambda1, lambda2 and eCprime
#' @keywords internal
gfl_get_value = function(valuehat, weight, trails, lambda1, lambda2, eCprime,
                         alpha=0.2, inflate=2, tol.beta=1e-6, maxsteps=100000) {
  #assume lambda1=0 and compute the fused lasso solution, centered around eCprime
  z=rep(0,tail(trails$breakpoints,n=1))
  u=rep(0,tail(trails$breakpoints,n=1))
  value = csnorm:::weighted_graphfl(valuehat, weight, trails$ntrails, trails$trails,
                         trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.beta, z, u) - eCprime
  #now soft-threshold the shifted value around eCprime
  value=sign(value)*pmax(abs(value)-lambda1, 0)
  return(value)
}

#' compute sparse fused lasso degrees of freedom
#' @keywords internal
get_gfl_degrees_of_freedom = function(mat, trails, tol.value=1e-3) {
  cl = csnorm:::build_patch_graph(mat, trails, tol.value=tol.value)$components
  mat[,patchno:=cl$membership]
  dof = mat[abs(value)>tol.value,uniqueN(patchno)] #sparse fused lasso
  stopifnot(dof<=cl$no)
  return(dof)
}

#' compute BIC for a given value of of lambda1, lambda2 and eCprime
#' @keywords internal
gfl_BIC = function(matg, trails, lambda1, lambda2, eCprime, tol.value=1e-3) {
  #get value with lambda1 set to zero to avoid round-off errors in degrees of freedom
  submat = matg[,.(name,bin1,bin2,valuehat,weight,ncounts)]
  submat[,value:=csnorm:::gfl_get_value(valuehat, weight, trails,
                                        lambda1=0, lambda2=lambda2, eCprime=eCprime)]
  #get the number of patches and deduce degrees of freedom
  cl = csnorm:::build_patch_graph(submat, trails, tol.value=tol.value)$components
  submat[,patchno:=cl$membership]
  dof = submat[abs(value)+tol.value>lambda1,uniqueN(patchno)] #sparse fused lasso
  stopifnot(dof<=cl$no)
  #now soft-threshold the value around eCprime
  submat[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
  #compute BIC
  BIC = submat[,sum(weight*((valuehat-(value+eCprime))^2))+log(.N)*dof]
  #compute mallow's Cp
  #Cp = submat[,sum(weight*((valuehat-(value+eCprime))^2 - 1))]+2*dof
  return(BIC)
}

#' cross-validate lambda1 and eCprime
#' 
#' @keywords internal
optimize_lambda1_eCprime = function(matg, trails, tol=1e-3, lambda2=0, eC.num=100) {
  #compute values for lambda1=0 and eCprime=0
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, lambda2, 0)]
  matg[,eCprime:=NULL]
  matg[,value.ori:=value]
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, trails, tol.value=tol)$components
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1]),by=patchno][,.N,keyby=value]
  stopifnot(patches[,sum(N)]==cl$no)
  #
  dt = foreach (lambda1=c(0,patches[,abs(value)],patches[,max(abs(value))]+tol),.combine=rbind) %:% 
    foreach (eCprime=seq(patches[,min(value)]-tol,patches[,max(value)]+tol,l=eC.num), .combine=rbind) %dopar% {
    matg[,value:=value.ori-eCprime]
    dof = matg[abs(value)>lambda1,uniqueN(patchno)] #sparse fused lasso
    stopifnot(dof<=cl$no)
    #now soft-threshold the value around eCprime
    matg[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
    #compute BIC
    BIC = matg[,sum(weight*((valuehat-(value+eCprime))^2))+log(.N)*dof]
    data.table(eCprime=eCprime,lambda1=lambda1,BIC=BIC,dof=dof)
  }
  #ggplot(dt)+geom_point(aes(eCprime,lambda1,fill=dof))
  #ggplot(dt[BIC<7000])+geom_line(aes(eCprime,BIC))+facet_wrap(~ lambda1)
  #ggplot(dt[BIC<7000])+geom_line(aes(lambda1,BIC,colour=BIC))+facet_wrap(~ eCprime)
  #ggplot(dt[lambda1<0.1])+geom_line(aes(lambda1,dof))
  #ggplot(dt)+geom_line(aes(lambda1,BIC))
  #ggplot(dt[lambda1<0.2])+geom_point(aes(lambda1,BIC))
  lambda1=dt[BIC==min(BIC),lambda1]
  eCprime=dt[BIC==min(BIC),eCprime]
  return(list(lambda1=lambda1,eCprime=eCprime))
}

#' cross-validate lambda1 and set eCprime at minimum
#' 
#' @keywords internal
optimize_lambda1_only = function(matg, trails, tol=1e-3, lambda2=0) {
  #set eCprime to lower bound when lambda1=0
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, lambda2, 0)]
  eCprime = matg[,min(value)-tol]
  matg[,value.ori:=value]
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, trails, tol.value=tol)$components
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1]),by=patchno][,.N,keyby=value]
  stopifnot(patches[,sum(N)]==cl$no)
  #
  dt = foreach (lambda1=c(0,patches[,value]-eCprime,patches[,max(value)]-eCprime+tol),.combine=rbind) %do% {
      matg[,value:=value.ori-eCprime]
      dof = matg[abs(value)>lambda1,uniqueN(patchno)] #sparse fused lasso
      stopifnot(dof<=cl$no)
      #now soft-threshold the value around eCprime
      matg[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
      #compute BIC
      BIC = matg[,sum(weight*((valuehat-(value+eCprime))^2))+log(.N)*dof]
      data.table(lambda1=lambda1,BIC=BIC,dof=dof)
    }
  #ggplot(dt)+geom_line(aes(lambda1,dof))
  #ggplot(dt[lambda1<0.1])+geom_line(aes(lambda1,dof))
  #ggplot(dt)+geom_line(aes(lambda1,BIC))
  #ggplot(dt[lambda1<0.2])+geom_point(aes(lambda1,BIC))
  lambda1=dt[BIC==min(BIC),lambda1][1]
  matg[,c("value","value.ori","patchno"):=NULL]
  return(list(lambda1=lambda1,eCprime=eCprime))
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2 = function(matg, trails, tol=1e-3, lambda1=0, eCprime=0, lambda2.last=1) {
  obj = function(x){csnorm:::gfl_BIC(matg, trails, lambda1=lambda1, lambda2=10^(x), eCprime=eCprime)}
  minlambda=tol/2
  maxlambda=10*lambda2.last
  #save(matg,trails,tol,lambda1,eCprime,file="debug_lambda2.RData")
  #cat("*** maxlambda ",maxlambda," (range=",matg[,max(value)-min(value)],")\n")
  #shrink maximum lambda, in case initial guess is too big
  repeat {
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, maxlambda, eCprime)]
    #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value)))
    if (matg[,max(value)-min(value)] > tol) {
      maxlambda = maxlambda*2
      break
    }
    maxlambda = maxlambda/2
    #dof = csnorm:::get_gfl_degrees_of_freedom(matg, trails)
    #cat("shrink maxlambda ",maxlambda," (range=",matg[,max(value)-min(value)],", dof=",dof,")\n")
  }
  #cat("optimize lambda2 between ",minlambda," and ",maxlambda,"\n")
  #dt = foreach (lam=seq(log10(minlambda),log10(maxlambda),l=100),.combine=rbind) %dopar% data.table(x=lam,y=obj(lam))
  #ggplot(dt)+geom_line(aes(10^(x),y))+scale_x_log10()
  op=optimize(obj, c(log10(minlambda),log10(maxlambda)), tol=tol)
  lambda2=10^op$minimum
  if (lambda2==minlambda | lambda2==maxlambda) cat("   Warning: lambda2 hit boundary.")
  if (lambda2<tol) lambda2=0
  return(lambda2)
}

#' run fused lasso on one dataset contained in matg, fusing 'value'
#' 
#' finds optimal lambda1, lambda2 and eC using BIC
#' @keywords internal
csnorm_fused_lasso_signal = function(matg, trails, tol=1e-3, verbose=T, ncores=ncores) {
  groupname=matg[,as.character(name[1])]
  eCprime=matg[,min(valuehat)]
  lambda1=0
  lambda2=1
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+scale_fill_gradient2())
  lambda2 = csnorm:::optimize_lambda2(matg, trails, tol=tol, lambda1=lambda1, eCprime=eCprime,
                                      lambda2.last=lambda2)
  #get best lambda1 and set eCprime to lower bound
  vals = csnorm:::optimize_lambda1_only(matg, trails, tol=tol, lambda2=lambda2)
  lambda1=vals$lambda1
  eCprime=vals$eCprime
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, lambda2, eCprime)]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  return(data.table(name=groupname,lambda1=lambda1,lambda2=lambda2,eCprime=eCprime))
}

#' run fused lasso on one dataset contained in matg, fusing 'value'
#' 
#' finds optimal lambda1, lambda2 and eC using BIC
#' @keywords internal
csnorm_fused_lasso_differential = function(matg, params, trails, tol=1e-3, verbose=T, ncores=ncores) {
  groupname=matg[,as.character(name[1])]
  if (is.null(params)) {
    eCprime=0
    lambda1=0
  } else {
    lambda1=params[name==groupname,lambda1]
    eCprime=params[name==groupname,eCprime]
  }
  matg[,value:=0]
  lambda2=1
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+scale_fill_gradient2())
  lambda2 = csnorm:::optimize_lambda2(matg, trails, tol=tol, lambda1=lambda1, eCprime=eCprime,
                                      lambda2.last=lambda2)
  #find best lambda1 and eCprime
  vals = csnorm:::optimize_lambda1_eCprime(matg, trails, tol=tol, lambda2=lambda2)
  lambda1=vals$lambda1
  eCprime=vals$eCprime
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, lambda2, eCprime)]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  return(data.table(name=groupname,lambda1=lambda1,lambda2=lambda2,eCprime=eCprime))
}

#' produce fused lasso solution at the given set of parameters
#' 
#' @keywords internal
get_lasso_coefs = function(matg, p, trails) {
  matg[,c("lambda1","lambda2","eCprime"):=list(p$lambda1, p$lambda2, p$eCprime)]
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, p$lambda1, p$lambda2, p$eCprime)]
  matg
}

#' produce fused lasso solution at the given set of parameters
#' 
#' @keywords internal
prepare_binless_decay = function(csg, mat) {
  #add new signal values to cts
  cts=merge(csg@cts[,.(name,bin1,bin2,dbin,cat,count,lmu.base,log_decay,eC,weight)],
            mat[,.(name,bin1,bin2,phi,eCprime)], by=c("name","bin1","bin2"), all.x=T)
  cts[,c("mu","kappaij"):=list(exp(eC+eCprime+lmu.base+log_decay+phi),eC+eCprime+log_decay)]
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/csg@par$alpha))]
  dbins=csg@par$dbins
  csd = cts[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
               kappahat=weighted.mean(z+kappaij, weight/var),
               std=1/sqrt(sum(weight/var)), weight=sum(weight)/2), keyby=c("name", "dbin")] #each count appears twice
  return(csd)
}

#' compute input to fused lasso
#' @keywords internal
csnorm_compute_raw_signal = function(csg, mat) {
  cts = copy(csg@cts)
  mat=mat[,.(name,bin1,bin2,phi)]
  cts = mat[cts,,on=c("name","bin1","bin2")]
  cts[,c("z","var"):=list(count/exp(phi+eC+lmu.base+log_decay)-1,
                          (1/exp(phi+eC+lmu.base+log_decay)+1/csg@par$alpha))]
  mat = cts[,.(phihat=weighted.mean(z+phi, weight/var),
               phihat.var=1/sum(weight/var),
               ncounts=sum(ifelse(count>0,weight,0))),keyby=c("name","bin1","bin2")][mat]
  mat[is.na(phihat),c("phihat","phihat.var","ncounts"):=list(1,Inf,0)] #bins with no detectable counts
  mat[,c("valuehat","weight"):=list(phihat,1/phihat.var)]
  setkey(mat,name,bin1,bin2)
  return(mat)
}

#' compute input to fused lasso
#' @keywords internal
csnorm_compute_raw_differential = function(csg, mat, ref) {
  ctsref = foreach(n=csg@cts[name!=ref,unique(name)],.combine=rbind) %do%
    csg@cts[name==ref,.(name=n,bin1,bin2,count,mu,weight)]
  cts=csg@cts[name!=ref]
  #
  if (is.null(mat)) {
    mat=cts[,CJ(name=name,bin1=ordered(levels(bin1)),bin2=ordered(levels(bin2)),sorted=T,unique=T)][bin2>=bin1]
    mat[,c("phi.ref","delta","diffsig"):=list(0,0,1)]
  } else {
    mat=mat[,.(name,bin1,bin2,phi.ref,delta,diffsig)]
  }
  ctsref=mat[ctsref]
  ctsref[,c("z","var"):=list(count/(exp(phi.ref)*mu)-1,(1/(exp(phi.ref)*mu)+1/csg@par$alpha))]
  mat.ref=ctsref[,.(phihat.ref=weighted.mean(z+phi.ref, weight/var),
                    sigmasq.ref=1/sum(weight/var),
                    ncounts.ref=sum(ifelse(count>0,weight,0))),keyby=c("name","bin1","bin2")][mat[,.(name,bin1,bin2)]]
  mat.ref[is.na(phihat.ref),c("phihat.ref","sigmasq.ref","ncounts.ref"):=list(1,Inf,0)] #bins with no detectable counts
  #
  cts=mat[cts]
  cts[,c("z","var"):=list(count/(exp(phi.ref+delta)*mu)-1,
                          (1/(exp(phi.ref+delta)*mu)+1/csg@par$alpha))]
  mat=cts[,.(phihat=weighted.mean(z+phi.ref+delta, weight/var),
             sigmasq=1/sum(weight/var),
             ncounts=sum(ifelse(count>0,weight,0))),keyby=c("name","bin1","bin2")][mat]
  mat[is.na(phihat),c("phihat","sigmasq","ncounts"):=list(1,Inf,0)] #bins with no detectable counts
  stopifnot(mat[,.N]==mat.ref[,.N])
  mat=merge(mat,mat.ref)
  mat[,c("deltahat","deltahat.var","ncounts"):=list(phihat-phihat.ref,sigmasq+sigmasq.ref,ncounts+ncounts.ref)]
  mat[,c("valuehat","weight","ncounts.ref"):=list(deltahat,1/deltahat.var,NULL)]
  setkey(mat,name,bin1,bin2)
  return(mat)
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
detect_binless_interactions = function(cs, resolution, group, fit.decay=F, ncores=1, niter=10, tol=1e-3,
                                       verbose=T, max_perf_iteration=1000, convergence_epsilon=1e-5){
  if (verbose==T) cat("Binless interaction detection with resolution=",resolution," and group=",group,"\n")
  if (group!="all" & fit.decay!=F) stop("Cannot yet fit decay with grouped data!")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csg, type="binteractions", threshold=-1, ref="expected", raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  #
  ### build matrix
  #create an empty matrix containing all cells, even those with no cut-site intersection
  bins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
  bins=unique(cut(c(bins,head(bins,n=-1)+resolution/2), bins, ordered_result=T, right=F, include.lowest=T,dig.lab=12))
  mat=CJ(name=csg@cts[,unique(name)],bin1=bins,bin2=bins,sorted=F,unique=F)[bin2>=bin1]
  mat[,phi:=0]
  ### build optimization trails
  if (verbose==T) cat("  Build triangle grid trails\n")
  trails = csnorm:::gfl_compute_trails(mat[,nlevels(bin1)])
  stopifnot(all(mat[,.N,by=name]$N==mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  stopifnot(all(length(V(trails$graph))==mat[,.N,by=name]$N))
  ### main loop
  params=NULL
  registerDoParallel(cores=ncores)
  for (step in 1:niter) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_signal(csg, mat)
    #
    #perform fused lasso on signal
    if (verbose==T) cat("  Fused lasso\n")
    groupnames=mat[,as.character(unique(name))]
    params = foreach(g=groupnames, .combine=rbind) %dopar%
      csnorm:::csnorm_fused_lasso_signal(mat[name==g], trails, tol=tol, ncores=ncores, verbose=verbose)
    #display param info
    if (verbose==T)
      for (i in 1:params[,.N])
        cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
            " eC'=",params[i,eCprime],"\n")
    #compute matrix at new params
    save(mat,params,file=paste0("mat_step_",step,".RData"))
    mat = foreach (g=groupnames, .combine=rbind) %dopar%
      csnorm:::get_lasso_coefs(mat[name==g],params[name==g], trails)
    #convert back value to the actual signal
    mat[,phi.old:=phi]
    mat[,phi:=value]
    cts=csg@cts
    if ("phi" %in% names(cts)) cts[,phi:=NULL]
    csg@cts=merge(cts, mat[,.(name,bin1,bin2,phi)], by=c("name","bin1","bin2"), all.x=T)
    #
    p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)
    ggsave(p,filename = paste0("sig_step_",step,"_value.png"), width=10, height=8)
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("sig_step_",step,"_weight.png"), width=10, height=8)
    #
    #check convergence
    if(mat[,all(abs(phi-phi.old)<tol)]) break
    #
    if (fit.decay!=F) {
      if (verbose==T) cat("  Fit decay\n")
      #compute diagonal decay
      csd = csnorm:::prepare_binless_decay(csg, mat)
      op = csnorm:::csnorm_gauss_decay_optimize(csd, csg@par$design, csg@par$Kdiag, csg@par$lambda_diag,
                                                'perf', max_perf_iteration=max_perf_iteration,
                                                convergence_epsilon=convergence_epsilon)
      p=ggplot(op$par$decay)+geom_point(aes(distance,kappahat))+geom_line(aes(distance,kappa))+
        scale_x_log10()+facet_wrap(~name)
      ggsave(p,filename = paste0("sig_step_",step,"_decay.png"), width=10, height=8)
      #report new parameters
      cts=csg@cts
      cts[,c("eC","log_decay"):=NULL]
      cts = merge(cbind(csg@par$design[,.(name)],eC=op$par$eC), cts, by="name", all.x=F,all.y=T)
      csg@cts = merge(op$par$decay[,.(name,dbin,log_decay)], cts, by=c("name","dbin"))
    }
  }
  #
  if (verbose==T) cat(" Detect patches\n")
  mat = csnorm:::detect_binless_patches(mat, trails, tol.value=tol)
  #
  ### store interaction
  csi=new("CSinter", mat=mat, type="binteractions", threshold=-1, ref="expected")
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
detect_binless_differences = function(cs, resolution, group, ref, niter=10, tol=1e-3, ncores=1, verbose=T){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=",group," and ref=",ref,"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="bdifferences", threshold=-1, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat=NULL
  trails=NULL
  params=NULL
  registerDoParallel(cores=ncores)
  for (step in 1:niter) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_differential(csg, mat, ref)
    #
    if (is.null(trails)) {
      if (verbose==T) cat("  Build triangle grid trails\n")
      trails = csnorm:::gfl_compute_trails(mat[,nlevels(bin1)])
    }
    #
    #perform fused lasso on signal
    if (verbose==T) cat("  Fused lasso\n")
    groupnames=mat[,as.character(unique(name))]
    params = foreach(g=groupnames, .combine=rbind) %dopar%
      csnorm:::csnorm_fused_lasso_differential(mat[name==g], params, trails, tol=tol, ncores=ncores, verbose=verbose)
    #display param info
    if (verbose==T)
      for (i in 1:params[,.N])
        cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
            " eC'=",params[i,eCprime],"\n")
    #compute matrix at new params
    #save(mat,params,file=paste0("dmat_step_",step,".RData"))
    mat = foreach (g=groupnames, .combine=rbind) %dopar%
      csnorm:::get_lasso_coefs(mat[name==g],params[name==g], trails)
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("diff_step_",step,"_value.png"), width=10, height=8)
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("diff_step_",step,"_weight.png"), width=10, height=8)
    #
    #convert back value to the actual signal
    mat[,delta:=value]
    mat[,diffsig.old:=diffsig]
    mat[,diffsig:=exp(delta)]
    mat[,phi.ref:=(phihat.ref/sigmasq.ref + (phihat-delta)/sigmasq)/(1/sigmasq.ref+1/sigmasq)]
    mat=mat[,.(name,bin1,bin2,phi.ref,delta,diffsig,diffsig.old,value,weight)]
    if(mat[,all(abs(diffsig-diffsig.old)<tol)]) break
  }
  if (verbose==T) cat(" Detect patches\n")
  mat = csnorm:::detect_binless_patches(mat, trails, tol.value=tol)
  csi=new("CSinter", mat=mat, type="bdifferences", threshold=-1, ref=ref)
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
