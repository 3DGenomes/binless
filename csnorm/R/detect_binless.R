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

#' compute trail information for the gfl package, for a triangle trid with nrow rows
#' @keywords internal
gfl_compute_trails = function(nrow) {
  chain = gfl_triangle_grid_chain(nrow)
  trails = gfl_chains_to_trails(chain)
  stopifnot(uniqueN(trails$trails)==nrow*(nrow+1)/2)
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
  #now soft-threshold manually around eCprime
  value=sign(value)*pmax(abs(value)-lambda1, 0)
  return(value)
}

#' compute BIC for a given value of of lambda1, lambda2 and eCprime
#' @keywords internal
gfl_BIC = function(matg, trails, lambda1, lambda2, eCprime, tol.value=1e-3) {
  #compute coefficients on each model
  submat = matg[,.(name,bin1,bin2,valuehat,weight,
                   value=csnorm:::gfl_get_value(valuehat, weight, trails,
                                                lambda1=lambda1, lambda2=lambda2, eCprime=eCprime))]
  #get the number of patches and degrees of freedom
  submat = csnorm:::detect_binless_patches(submat, tol.value=tol.value)
  dof = submat[abs(value)>tol.value,uniqueN(patchno)] #sparse fused lasso
  #compute BIC
  BIC = submat[,sum(weight*((valuehat-(value+eCprime))^2))+log(.N)*dof]
  #compute mallow's Cp
  #Cp = submat[,sum(weight*((valuehat-(value+eCprime))^2 - 1))]+2*dof
  return(BIC)
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

#' give a patch ID to each patch in a binless matrix, and report local extrema
#' @keywords internal
detect_binless_patches = function(mat, tol.value=1e-3) {
  if (mat[,uniqueN(name)]>1)
    return(foreach (n=mat[,unique(name)],.combine=rbind) %do% csnorm:::detect_binless_patches(mat[name==n]))
  #build bin graph
  g = csnorm:::compute_2d_connectivity_graph(mat[,as.integer(max(bin1))])
  stopifnot(length(V(g))==mat[,.N])
  V(g)$value = mat[,value]
  V(g)$weight = 1
  #plot(g, vertex.color=factor(V(g)$value), vertex.label=V(g)$count,
  #     layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #
  #remove edges that connect different values
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
  #
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

#' compute input to fused lasso
#' @keywords internal
csnorm_compute_raw_signal = function(csg, mat) {
  cts = csg@cts
  if (is.null(mat)) mat=cts[,.(phi=0,signal=1,eCprime=0),by=c("name","bin1","bin2")]
  cts = mat[,.(name,bin1,bin2,phi,signal,eCprime)][cts,,on=c("name","bin1","bin2")]
  cts[,c("z","var"):=list(count/(exp(phi+eCprime)*mu)-1,(1/(exp(phi+eCprime)*mu)+1/csg@dispersion))]
  mat = cts[,.(phihat=weighted.mean(z+phi+eCprime, weight/var),
               phihat.var=1/sum(weight/var)),by=c("name","bin1","bin2","phi","signal")]
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
  if (is.null(mat)) mat=cts[,.(phi.ref=0,delta=0,diffsig=1),by=c("name","bin1","bin2")]
  ctsref=mat[ctsref]
  ctsref[,c("z","var"):=list(count/(exp(phi.ref)*mu)-1,(1/(exp(phi.ref)*mu)+1/csg@dispersion))]
  mat=mat[ctsref[,.(phihat.ref=weighted.mean(z+phi.ref, weight/var),
                    sigmasq.ref=1/sum(weight/var)),keyby=c("name","bin1","bin2")]]
  #
  cts=mat[cts]
  cts[,c("z","var"):=list(count/(exp(phi.ref+delta)*mu)-1,
                          (1/(exp(phi.ref+delta)*mu)+1/csg@dispersion))]
  mat=mat[cts[,.(phihat=weighted.mean(z+phi.ref+delta, weight/var),
                 sigmasq=1/sum(weight/var)),by=c("name","bin1","bin2")]]
  mat[,c("deltahat","deltahat.var"):=list(phihat-phihat.ref,sigmasq+sigmasq.ref)]
  mat[,c("valuehat","weight"):=list(deltahat,1/deltahat.var)]
  setkey(mat,name,bin1,bin2)
  return(mat)
}

#' cross-validate lambda1
#' @keywords internal
optimize_lambda1 = function(matg, trails, tol=1e-3, lambda2=0, eCprime=0, enforce.positivity=T) {
  obj = function(x){csnorm:::gfl_BIC(matg, trails, lambda1=x, lambda2=lambda2, eCprime=eCprime)}
  if (enforce.positivity==T) { #we force all signal values to be positive
    minlambda=matg[,abs(min(pmin(value,0)))] #in this case minlambda <= maxlambda always
    maxlambda=matg[,abs(max(value))]
    stopifnot(minlambda<=maxlambda)
  } else { #we don't
    minlambda=0
    maxlambda=matg[,max(abs(value))]
  }
  #dt = foreach (lam=seq(minlambda,maxlambda,l=100),.combine=rbind) %dopar% data.table(x=lam,y=obj(lam))
  #ggplot(dt)+geom_line(aes(x,y))
  if (maxlambda>0) {
    lambda1=optimize(obj, c(minlambda,maxlambda), tol=tol)$minimum
    if (abs(lambda1)<tol) lambda1=0
  } else {
    lambda1=0
  }
  return(lambda1)
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2 = function(matg, trails, tol=1e-3, lambda1=0, eCprime=0) {
  obj = function(x){csnorm:::gfl_BIC(matg, trails, lambda1=lambda1, lambda2=x, eCprime=eCprime)}
  minlambda=0
  maxlambda=matg[,max(valuehat)-min(valuehat)]
  #first shrink maximum lambda, in case initial guess is too big
  repeat {
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, maxlambda, eCprime)]
    #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value)))
    if (matg[,max(value)-min(value)] > tol) break
    maxlambda = maxlambda/2
  }
  #now expand it until we fuse everyone
  repeat {
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, maxlambda, eCprime)]
    #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value)))
    if (matg[,max(value)-min(value)] <= tol) break
    maxlambda = maxlambda*2
  }
  #dt = foreach (lam=seq(minlambda,maxlambda,l=100),.combine=rbind) %dopar% data.table(x=lam,y=obj(lam))
  #ggplot(dt)+geom_line(aes(x,y))
  op=optimize(obj, c(minlambda,maxlambda), tol=tol)
  return(op$minimum)
}

#' run fused lasso on each dataset contained in mat, fusing 'value'
#' 
#' finds optimal lambda1, lambda2 and eC through k-fold cv.
#' @keywords internal
csnorm_fused_lasso = function(mat, trails, tol=1e-3, niter=10, verbose=T, enforce.positivity=T, ncores=ncores) {
  groupnames=mat[,as.character(unique(name))]
  params = foreach(g=groupnames, .combine=rbind) %dopar% {
    matg=mat[name==g]
    #lambda2
    if (enforce.positivity==T) {
      eCprime=matg[,min(valuehat)]
      lambda1=matg[,mad(valuehat)]
    } else {
      eCprime=0
      lambda1=0
    }
    matg[,value:=0]
    #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+scale_fill_gradient2())
    for (i in 1:niter) {
      #lambda2=10
      #beta=csnorm:::weighted_graphfl(matg[,valuehat], matg[,weight], trails$ntrails, trails$trails, trails$breakpoints,
      #                               lambda2, 0.2, 2, 100000, 1e-3, rep(0,tail(trails$breakpoints,n=1)), rep(0,tail(trails$breakpoints,n=1)))
      #write.table(matg[,valuehat],file="tmp/gfl-master/example/mat_y.csv", quote=F, row.names=F, col.names=F)
      #write.table(matg[,weight],file="tmp/gfl-master/example/mat_weight.csv", quote=F, row.names=F, col.names=F)
      #matg[,value:=beta]
      #matg[,value.py:=fread("tmp/gfl-master/example/beta.csv")$V1]
      #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+geom_raster(aes(bin2,bin1,fill=value.py))
      #ggplot(matg)+geom_point(aes(value,value.py))
      #
      matg[,value.old:=value]
      lambda2 = csnorm:::optimize_lambda2(matg, trails, tol=tol, lambda1=lambda1, eCprime=eCprime)
      matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, lambda2, eCprime)]
      #lambda1
      lambda1 = csnorm:::optimize_lambda1(matg, trails, tol=tol, lambda2=lambda2, eCprime=eCprime,
                                          enforce.positivity=enforce.positivity)
      matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, lambda2, eCprime)]
      #eCprime
      if (enforce.positivity==T) {
        eCprime = matg[,min(value+eCprime)]
        matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1, lambda2, eCprime)]
      }
      #
      cat("   iteration ",i," : lambda1=",lambda1," lambda2=",lambda2," eC'=",eCprime,"\n")
      #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
      if (matg[,all(abs(value-value.old)<tol)]) break
    }
    data.table(name=g,lambda1=lambda1,lambda2=lambda2,eCprime=eCprime)
  }
  #
  if (verbose==T)
    for (i in 1:params[,.N])
      cat("   ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
          " eC'=",params[i,eCprime],"\n")
  #
  mat = foreach (g=groupnames, .combine=rbind) %do% {
    p=params[name==g]
    matg=mat[name==g]
    matg[,c("lambda1","lambda2","eCprime"):=list(p$lambda1, p$lambda2, p$eCprime)]
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, p$lambda1, p$lambda2, p$eCprime)]
    matg
  }
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2()+facet_wrap(~name)
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
detect_binless_interactions = function(cs, resolution, group, ncores=1,
                                       niter.outer=10, niter.inner=10, tol=1e-3, verbose=T){
  if (verbose==T) cat("Binless interaction detection with resolution=",resolution," and group=",group,"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csg, type="binteractions", threshold=-1, ref="expected", raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  #
  ### build matrix
  mat=NULL
  trails=NULL
  for (step in 1:niter.outer) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_signal(csg, mat)
    #
    if (is.null(trails)) {
      if (verbose==T) cat("  Build triangle grid trails\n")
      trails = csnorm:::gfl_compute_trails(mat[,nlevels(bin1)])
      stopifnot(all(mat[,.N,by=name]$N==mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
    }
    #perform fused lasso on signal to noise
    if (verbose==T) cat("  Fused lasso\n")
    mat = csnorm:::csnorm_fused_lasso(mat, trails, tol=tol,
                                      niter=niter.inner, enforce.positivity=T, ncores=ncores)
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=value))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")
    #
    #convert back value to the actual signal
    mat[,phi:=value]
    mat[,signal.old:=signal]
    mat[,signal:=exp(phi)]
    mat=mat[,.(name,bin1,bin2,phi,eCprime,signal,signal.old,value,valuehat,weight)]
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=log10(signal)))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")+geom_raster(aes(bin2,bin1,fill=log10(signal.old)))
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=log10(weight)))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")+geom_raster(aes(bin2,bin1,fill=log10(signal.old)))
    if(mat[,all(abs(signal-signal.old)<tol)]) break
  }
  if (verbose==T) cat(" Detect patches\n")
  mat = csnorm:::detect_binless_patches(mat)
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
detect_binless_differences = function(cs, resolution, group, ref, niter.outer=10, niter.inner=10,
                                      tol=1e-3, ncores=1, verbose=T){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=",group," and ref=",ref,"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="bdifferences", threshold=-1, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat=NULL
  trails=NULL
  for (step in 1:niter.outer) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_differential(csg, mat, ref)
    #
    if (is.null(trails)) {
      if (verbose==T) cat("  Build triangle grid trails\n")
      trails = csnorm:::gfl_compute_trails(mat[,nlevels(bin1)])
    }
    #
    #perform fused lasso on signal to noise
    if (verbose==T) cat("  Fused lasso\n")
    mat = csnorm:::csnorm_fused_lasso(mat, trails, tol=tol,
                                      niter=niter.inner, enforce.positivity=F, ncores=ncores)
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=value))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")
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
  mat = csnorm:::detect_binless_patches(mat)
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
plot_binless_matrix = function(mat, minima=F) {
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
      guides(fill=F)+scale_fill_gradient2()+facet_wrap(~name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="blue",fill=NA,data=b)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="red",fill=NA,data=a)
  } else {
    p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=-value))+
      geom_raster(aes(begin2,begin1,fill=-value))+
      guides(fill=F)+scale_fill_gradient2()+facet_wrap(~name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="black",fill=NA,data=a)
  }
  print(p)
}
