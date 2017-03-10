#' @include csnorm.R
NULL

#' connectivity on a triangle
#' @keywords internal
flsa_compute_connectivity = function(nbins, start=0) {
  stopifnot(start>=0,nbins>=2)
  if (nbins==2) {
    ret=sapply(list(start+1,c(start,start+2),start+1),as.integer)
    names(ret) = c(start, start+1, start+2)
  } else {
    upper.row = c(list(start+1),
                  sapply((start+1):(start+nbins-2),function(x){c(x-1,x+1,x+nbins-1)},simplify=F),
                  list(c(start+nbins-2,start+2*(nbins-1))))
    names(upper.row) = start:(start+nbins-1)
    lower.tri = flsa_compute_connectivity(nbins-1,start=start+nbins)
    for (i in 1:(nbins-1))
      lower.tri[[i]] = c(lower.tri[[i]], start+i)
    ret = sapply(c(upper.row,lower.tri),as.integer)
  }
  class(ret) = "connListObj"
  return(ret)
}


#' compute cv error for a given value of of lambda1 and lambda2
#' @keywords internal
flsa_get_value = function(flsa.op, lambda1, lambda2, eCsd) {
  #assume lambda1=0
  value=flsaGetSolution(flsa.op, lambda1=0, lambda2=lambda2)[1,]-eCsd
  #now soft-threshold manually around eCsd
  value=sign(value)*pmax(abs(value)-lambda1, 0)
  return(value)
}

#' compute Mallow's Cp for a given value of of lambda1, lambda2 and eCsd
#' @keywords internal
flsa_Cp = function(mat, res.ref, lambda1, lambda2, eCsd) {
  #compute coefficients on each model
  submat = mat[,.(name,bin1,bin2,ncounts,valuehat,
                  value=csnorm:::flsa_get_value(res.ref[[as.character(name[1])]],
                                                lambda1=lambda1, lambda2=lambda2, eCsd=eCsd))]
  #get the number of patches and degrees of freedom
  submat = csnorm:::detect_binless_patches(submat)
  dof = submat[value!=0,uniqueN(patchno)]
  #compute mallow's Cp
  Cp = submat[,sum(((valuehat-(value+eCsd))^2 - 1))]+2*dof
  return(Cp)
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
detect_binless_patches = function(mat) {
  if (mat[,uniqueN(name)]>1)
    return(foreach (n=mat[,unique(name)],.combine=rbind) %do% csnorm:::detect_binless_patches(mat[name==n]))
  #build bin graph
  g = compute_2d_connectivity_graph(mat[,as.integer(max(bin1))])
  stopifnot(length(V(g))==mat[,.N])
  V(g)$value = mat[,value]
  V(g)$weight = 1
  #plot(g, vertex.color=factor(V(g)$value), vertex.label=V(g)$count,
  #     layout=as.matrix(mat[,.(as.integer(bin1),as.integer(bin2))]))
  #
  #remove edges that connect different values
  values_by_edge = t(matrix(get.vertex.attribute(g,"value",index=get.edges(g,E(g))), ncol=2))
  to_rm = E(g)[abs(values_by_edge[1,]-values_by_edge[2,])>1e-8]
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
               phihat.sd=sqrt(1/sum(weight/var)),ncounts=sum(weight)),by=c("name","bin1","bin2","phi","signal")]
  mat[,valuehat:=phihat/phihat.sd]
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
  if (is.null(mat)) mat=cts[,.(phi.ref=0,delta=0,diffsig=1,ncounts=sum(weight)),by=c("name","bin1","bin2")]
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
  mat[,c("deltahat","deltahat.sd"):=list(phihat-phihat.ref,sqrt(sigmasq+sigmasq.ref))]
  mat[,valuehat:=deltahat/deltahat.sd]
  setkey(mat,name,bin1,bin2)
  return(mat)
}

#' cross-validate lambda1
#' @keywords internal
optimize_lambda1 = function(mat, res.ref, g, lambda2=0, eCsd=0, enforce.positivity=T) {
  obj = function(x){csnorm:::flsa_Cp(mat[name==g], res.ref, lambda1=x, lambda2=lambda2, eCsd=eCsd)}
  if (enforce.positivity==T) { #we force all signal values to be positive
    minlambda=mat[name==g,abs(min(value))] #in this case minlambda <= maxlambda always
    maxlambda=mat[name==g,abs(max(value))]
    stopifnot(minlambda<=maxlambda)
  } else { #we don't
    minlambda=0
    maxlambda=mat[name==g,max(abs(value))]
  }
  #ggplot(data.table(x=seq(minlambda,maxlambda,l=100))[,.(x,y=sapply(x,obj))])+geom_line(aes(x,y))
  if (maxlambda>0) {
    lambda1=optimize(obj, c(minlambda,maxlambda), tol=tol)$minimum
    if (abs(lambda1)<tol) lambda1=0
  } else {
    lambda1=0
  }
  data.table(name=g, lambda2=lambda2, eCsd=eCsd, lambda1=lambda1)
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2 = function(mat, res.ref, g, lambda1=0, eCsd=0) {
  obj = function(x){csnorm:::flsa_Cp(mat[name==g], res.ref, lambda1=lambda1, lambda2=x, eCsd=eCsd)}
  minlambda=0
  maxlambda=max(res.ref[[g]]$BeginLambda)
  #ggplot(data.table(x=seq(minlambda,maxlambda,l=100))[,.(x,y=sapply(x,obj))])+geom_line(aes(x,y))
  op=optimize(obj, c(minlambda,maxlambda), tol=tol)
  l2vals=sort(res.ref[[g]]$BeginLambda)
  i=findInterval(op$minimum, l2vals)
  if (i==length(l2vals)) {
    lambda2=l2vals[i]
  } else if (i==0) {
    lambda2=0
  } else if (obj(l2vals[i]) <= obj(l2vals[i+1])) {
    lambda2=l2vals[i]
  } else {
    lambda2=l2vals[i+1]
  }
  data.table(name=g, lambda2=lambda2, eCsd=eCsd, lambda1=lambda1)
}

#' run fused lasso on each dataset contained in mat, fusing 'value'
#' 
#' finds optimal lambda1, lambda2 and eC through k-fold cv.
#' @keywords internal
csnorm_fused_lasso = function(mat, tol=1e-3, niter=10, verbose=T, enforce.positivity=T, ncores=ncores) {
  if (verbose==T) cat("   compute fused lasso solution path\n")
  groupnames=mat[,as.character(unique(name))]
  cmat = csnorm:::flsa_compute_connectivity(mat[,nlevels(bin1)])
  registerDoParallel(cores=ncores)
  res.ref = foreach (g=groupnames, .final=function(x){setNames(x,groupnames)}) %dopar% {
    submat=mat[name==g]
    stopifnot(length(cmat)==submat[,.N]) #ibin indices have something wrong
    flsa(submat[,valuehat], connListObj=cmat, verbose=F)
  }
  #fail if any has failed (flsa bugs)
  if (any(sapply(res.ref, is.null))) stop("one flsa run failed, aborting")
  #
  #determine optimal parameters
  if (verbose==T) cat("   determine optimal parameters\n")
  params = foreach(g=groupnames, .combine=rbind) %dopar% {
    matg=mat[name==g]
    #lambda2
    lambda1=matg[,mad(valuehat)]
    if (enforce.positivity==T) {
      eCsd=matg[,min(valuehat)]
    } else {
      eCsd=0
    }
    matg[,value:=0]
    for (i in 1:niter) {
      matg[,value.old:=value]
      lambda2 = csnorm:::optimize_lambda2(matg, res.ref, g,
                                          lambda1=lambda1, eCsd=eCsd)[,lambda2]
      matg[,value:=csnorm:::flsa_get_value(res.ref[[g]], lambda1=lambda1, lambda2=lambda2, eCsd=eCsd)]
      #lambda1
      lambda1 = csnorm:::optimize_lambda1(matg, res.ref, g,
                                          lambda2=lambda2, eCsd=eCsd,
                                          enforce.positivity=enforce.positivity)[,lambda1]
      matg[,value:=csnorm:::flsa_get_value(res.ref[[g]], lambda1=lambda1, lambda2=lambda2, eCsd=eCsd)]
      #eCsd
      if (enforce.positivity==T) {
        eCsd = matg[,min(value+eCsd)]
        matg[,value:=csnorm:::flsa_get_value(res.ref[[g]], lambda1=lambda1, lambda2=lambda2, eCsd=eCsd)]
      }
      #
      cat("   iteration ",i," : lambda1=",lambda1," lambda2=",lambda2," eCsd'=",eCsd,"\n")
      #ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2()
      if (matg[,all(abs(value-value.old)<tol)]) break
    }
    data.table(name=g,lambda1=lambda1,lambda2=lambda2,eCsd=eCsd)
  }
  mat = foreach (g=groupnames, .combine=rbind) %do% {
    p=params[name==g]
    matg=mat[name==g]
    matg[,c("lambda1","lambda2","eCsd"):=list(p$lambda1, p$lambda2, p$eCsd)]
    matg[,value:=csnorm:::flsa_get_value(res.ref[[g]], lambda1=p$lambda1, lambda2=p$lambda2,eCsd=p$eCsd)]
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
#' @param niter number of IRLS iterations, and Cp iterations within
#'   
#' @return 
#' @export
#' 
#' @examples
detect_binless_interactions = function(cs, resolution, group, ncores=1, niter=10, tol=1e-3, verbose=T){
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
  for (step in 1:niter) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_signal(csg, mat)
    #
    #perform fused lasso on signal to noise
    if (verbose==T) cat("  Fused lasso\n")
    mat = csnorm:::csnorm_fused_lasso(mat, tol=tol, niter=niter, enforce.positivity=T, ncores=ncores)
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=value))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")
    #
    #convert back value to the actual signal
    mat[,c("phi","eCprime"):=list(value*phihat.sd,eCsd*phihat.sd)]
    mat[,signal.old:=signal]
    mat[,signal:=exp(phi)]
    mat=mat[,.(name,bin1,bin2,phi,eCprime,signal,signal.old,ncounts,value,valuehat)]
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=log10(signal)))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")+geom_raster(aes(bin2,bin1,fill=log10(signal.old)))
    if(mat[,all(abs(signal-signal.old)<tol)]) break
  }
  mat[,ncounts:=NULL]
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
detect_binless_differences = function(cs, resolution, group, ref, niter=10, tol=1e-3, ncores=1, verbose=T){
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                      " group=",group," and ref=",ref,"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="bdifferences", threshold=-1, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat=NULL
  for (step in 1:niter) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_differential(csg, mat, ref)
    #
    #perform fused lasso on signal to noise
    if (verbose==T) cat("  Fused lasso\n")
    mat = csnorm:::csnorm_fused_lasso(mat, tol=tol, niter=niter, enforce.positivity=F, ncores=ncores)
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=value))+facet_wrap(~name)+scale_fill_gradient(na.value = "black")
    #
    #convert back value to the actual signal
    mat[,delta:=value*deltahat.sd]
    mat[,diffsig.old:=diffsig]
    mat[,diffsig:=exp(delta)]
    mat[,phi.ref:=(phihat.ref/sigmasq.ref + (phihat-delta)/sigmasq)/(1/sigmasq.ref+1/sigmasq)]
    mat=mat[,.(name,bin1,bin2,phi.ref,delta,diffsig,diffsig.old,ncounts,value)]
    if(mat[,all(abs(diffsig-diffsig.old)<tol)]) break
  }
  mat[,ncounts:=NULL]
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
