#' @include csnorm.R
NULL

#' common calculation between interactions and differences
#' @keywords internal
estimate_binless_common_irls = function(cs, resolution, group, ref="expected", ncores=1, verbose=T) {
  # get zeros per bin
  if (verbose==T) cat("  Get zeros per bin\n")
  zeros = csnorm:::get_nzeros_binning(cs, resolution, ncores=ncores)
  # predict means
  if (verbose==T) cat("  Predict means\n")
  cts = csnorm:::csnorm_predict_binned_muhat_irls(cs, resolution, zeros)
  # group
  if (verbose==T) cat("  Group\n")
  if (group=="all") {
    names=cs@experiments[,unique(name)]
    groups=data.table(name=names,groupname=names)
  } else {
    groups=cs@experiments[,.(name,groupname=do.call(paste,mget(group))),by=group][,.(name,groupname)] #we already know groupname is unique
    groups[,groupname:=ordered(groupname)] #same class as name
  }
  setkey(groups,name)
  cts = groups[cts]
  cts[,name:=NULL]
  setnames(cts,"groupname","name")
  #
  if (ref != "expected") {
    if (!(ref %in% groups[,groupname]))
      stop("Reference group name not found! Valid names are: ",deparse(as.character(groups[,groupname])))
    if (groups[groupname!=ref,.N]==0)
      stop("There is no other group than ",ref, ", cannot compute differences!")
  }
  setkeyv(cts,c("name","bin1","bin2"))
  return(cts)
}

#' group counts to compute signal by cross-validated lasso
#' @keywords internal
estimate_binless_signal_irls = function(cs, resolution, group, mat, verbose=T) {
  cts=csnorm:::estimate_binless_common_irls(cs, resolution, group, ref="expected",
                                       ncores=ncores, verbose=verbose)
  if (verbose==T) cat("  Phi\n")
  init=cs@par
  if (is.null(mat)) mat=cts[,.(phi=0,signal=1),by=c("name","bin1","bin2")]
  cts=mat[,.(name,bin1,bin2,phi,signal)][cts,,on=c("name","bin1","bin2")]
  cts[,c("z","var"):=list(count/(exp(phi)*mu)-1,(1/(exp(phi)*mu)+1/init$alpha))]
  ret=cts[,.(phihat=weighted.mean(z+phi, weight/var),
             phihat.sd=sqrt(1/sum(weight/var)),ncounts=sum(weight)),by=c("name","bin1","bin2","phi","signal")]
  ret[,value:=phihat/phihat.sd]
  return(ret)
}

#' group counts to compute signal by cross-validated lasso
#' @keywords internal
estimate_binless_differential_irls = function(cs, resolution, group, mat, ref, verbose=T) {
  cts=csnorm:::estimate_binless_common_irls(cs, resolution, group, ref=ref,
                                            ncores=ncores, verbose=verbose)
  if (verbose==T) cat("  Delta\n")
  init=cs@par
  ctsref = foreach(n=cts[name!=ref,unique(name)],.combine=rbind) %do%
    cts[name==ref,.(name=n,bin1,bin2,count,mu,weight)]
  cts=cts[name!=ref]
  #
  if (is.null(mat)) mat=cts[,.(phi.ref=0,delta=0,diffsig=1,ncounts=sum(weight)),by=c("name","bin1","bin2")]
  ctsref=mat[ctsref]
  ctsref[,c("z","var"):=list(count/(exp(phi.ref)*mu)-1,(1/(exp(phi.ref)*mu)+1/init$alpha))]
  mat=mat[ctsref[,.(phihat.ref=weighted.mean(z+phi.ref, weight/var),
                    sigmasq.ref=1/sum(weight/var)),keyby=c("name","bin1","bin2")]]
  #
  cts=mat[cts]
  cts[,c("z","var"):=list(count/(exp(phi.ref+delta)*mu)-1,
                          (1/(exp(phi.ref+delta)*mu)+1/init$alpha))]
  mat=mat[cts[,.(phihat=weighted.mean(z+phi.ref+delta, weight/var),
                 sigmasq=1/sum(weight/var)),by=c("name","bin1","bin2")]]
  mat[,c("deltahat","deltahat.sd"):=list(phihat-phihat.ref,sqrt(sigmasq+sigmasq.ref))]
  mat[,value:=deltahat/deltahat.sd]
  return(mat)
}

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
flsa_cross_validate = function(mat, res.cv, cvgroups, lambda1, lambda2) {
  #compute coefficients on each model
  submat.all = mat[,.(name,bin1,bin2,ncounts, value.ref=value)]
  #compute cv coefs for each group, given lambda
  submat.cv = foreach(cv=cvgroups, .combine=rbind) %do% {
    ret = mat[,.(name,cv.group,bin1,bin2)]
    if (lambda1>0) {
      ret[,value:=flsaGetSolution(res.cv[[cv]][[as.character(name[1])]],
                                        lambda1=lambda1, lambda2=lambda2)[1,1,],by=name]
    } else {
      ret[,value:=flsaGetSolution(res.cv[[cv]][[as.character(name[1])]],
                                        lambda1=lambda1, lambda2=lambda2)[1,],by=name]
    }
    ret[cv.group==cv]
  }
  #ggplot(submat.cv)+geom_raster(aes(ibin1,ibin2,fill=value))+facet_grid(groupname~cv.group)
  #ggplot(submat.all)+geom_raster(aes(ibin1,ibin2,fill=value.ref))+facet_grid(~name)#+scale_fill_gradient2()
  ret=merge(submat.cv, submat.all, all.x=T, by=c("name","bin1","bin2"), sort=F)
  ret[,mse:=ncounts*(value-value.ref)^2]
  ret=ret[,.(mse=mean(mse),ncounts=.N),by=c("name","cv.group")]
  ret=ret[,.(mse=weighted.mean(mse,ncounts),mse.sd=sd(mse)),by=name]
  ret[,.(name,lambda1=lambda1,lambda2=lambda2,mse,mse.sd)]
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

#' run fused lasso on each dataset contained in mat, fusing 'value'
#' 
#' finds optimal lambda2 through k-fold cv. Uses lambda1=0
#' @keywords internal
csnorm_fused_lasso = function(mat, cv.fold=10, cv.gridsize=30, verbose=T) {
  if (verbose==T) cat("  estimation runs\n")
  groupnames=mat[,as.character(unique(name))]
  cmat = csnorm:::flsa_compute_connectivity(mat[,nlevels(bin1)])
  res.ref = foreach (g=groupnames, .final=function(x){setNames(x,groupnames)}) %dopar% {
    submat=mat[name==g]
    stopifnot(length(cmat)==submat[,.N]) #ibin indices have something wrong
    flsa(submat[,value], connListObj=cmat, verbose=F)
  }
  #fail if any has failed (flsa bugs)
  if (any(sapply(res.ref, is.null))) stop("one flsa run failed, aborting")
  #cross-validate lambda2 (lambda1=0 gives good results)
  mat[,cv.group:=as.character(((unclass(bin2)+unclass(bin1)*max(unclass(bin1)))%%cv.fold)+1)]
  cvgroups=as.character(1:cv.fold)
  #first, run lasso regressions
  if (verbose==T) cat("  cross-validation runs\n")
  res.cv = foreach (cv=cvgroups, .final=function(x){setNames(x,cvgroups)}) %:%
    foreach (g=groupnames, .final=function(x){setNames(x,groupnames)}) %dopar% {
      submat=copy(mat[name==g])
      submat[cv.group==cv,value:=0]
      stopifnot(length(cmat)==submat[,.N]) #ibin indices have something wrong
      flsa(submat[,value], connListObj=cmat, verbose=F)
    }
  #remove failed cv attempts (flsa bugs)
  cvgroups = foreach(cv=cvgroups, .combine=c) %do% {
    if (any(sapply(res.cv[[cv]], is.null)))
      cat("removing CV group ",cv," due to failure of flsa!\n")
    else
      cv
  }
  #then, try a few lambda values to define the interesting region
  if (verbose==T) cat("  coarse scan\n")
  cv = foreach (g=groupnames, .combine=rbind) %do% {
    minlambda=min(res.ref[[g]]$EndLambda[res.ref[[g]]$EndLambda>0])
    maxlambda=max(res.ref[[g]]$BeginLambda)
    l2vals=c(0,10^seq(log10(minlambda),log10(maxlambda),length.out=cv.gridsize))
    foreach (lambda2=l2vals, .combine=rbind) %dopar%
      csnorm:::flsa_cross_validate(mat[name==g], res.cv, cvgroups, lambda1=0, lambda2=lambda2)
  }
  #ggplot(cv)+geom_line(aes(lambda2,mse))+geom_errorbar(aes(lambda2,ymin=mse-mse.sd,ymax=mse+mse.sd))+
  #  facet_wrap(~name,scales="free")+scale_x_log10()
  #lower bound is lambda one step below minimum mse
  bounds=merge(cv,cv[,.SD[mse==min(mse),.(lambda2.max=lambda2)],by=name],by="name")
  lmin=bounds[lambda2<lambda2.max,.(l2min=max(lambda2)),by=name]
  #upper bound is lambda one step above mse
  lmax=bounds[,.SD[lambda2>lambda2.max,.(l2max=min(lambda2))],by=name]
  bounds=merge(lmin,lmax,by="name", all=T)
  bounds[,l2max:=pmin(l2max, sapply(as.character(name),FUN=function(g){max(res.ref[[g]]$BeginLambda)}), na.rm=T)]
  #
  #now, try all lambda values between these bounds
  if (verbose==T) cat("  fine scan\n")
  cv = foreach (g=groupnames, .combine=rbind) %do% {
    l2vals=sort(unique(res.ref[[g]]$BeginLambda))
    l2vals=bounds[name==g,l2vals[l2vals>=l2min & l2vals<=l2max]]
    foreach (lambda2=l2vals, .combine=rbind) %dopar%
      csnorm:::flsa_cross_validate(mat[name==g], res.cv, cvgroups, lambda1=0, lambda2=lambda2)
  }
  #ggplot(cv)+geom_line(aes(lambda2,mse))+#geom_errorbar(aes(lambda2,ymin=mse-mse.sd,ymax=mse+mse.sd))+
  #  facet_wrap(~name,scales="free")
  #determine optimal lambda value and compute coefficients
  if (verbose==T) cat("  report optimum\n")
  l2vals=cv[,.SD[mse==min(mse),.(lambda2)],by=name]
  for (g in groupnames)
    mat[name==g,value:=flsaGetSolution(res.ref[[g]], lambda1=0, lambda2=l2vals[name==g,lambda2])[1,]]
  return(mat)
}

#' Perform binless interaction detection using fused lasso
#'
#' @param cs 
#' @param ref 
#' @param resolution 
#' @param group 
#' @param ncores 
#' @param niter number of IRLS iterations
#' @param cv.fold perform x-fold cross-validation
#' @param cv.gridsize compute possible values of log lambda on a grid
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
csnorm_detect_binless = function(cs, resolution, group, ref="expected", niter=10, tol=1e-3, ncores=1,
                                 cv.fold=10, cv.gridsize=30, verbose=T) {
  if (verbose==T) cat("Binless detection with resolution=",resolution," group=",group," and ref=",ref,"\n")
  mat=NULL
  for (step in 1:niter) {
    if (verbose==T) cat("Main loop, step ",step,"\n")
    if (verbose==T) cat(" Estimate raw signal\n")
    if (ref=="expected") {
      mat = csnorm:::estimate_binless_signal_irls(cs, resolution, group, mat, verbose)
    } else {
      mat = csnorm:::estimate_binless_differential_irls(cs, resolution, group, mat, ref, verbose)
    }
    setkey(mat,name,bin1,bin2)
    #
    #perform fused lasso on signal to noise
    if (verbose==T) cat(" Fused lasso\n")
    mat = csnorm:::csnorm_fused_lasso(mat, cv.fold=cv.fold, cv.gridsize=cv.gridsize)
    #ggplot(mat)+geom_raster(aes(ibin1,ibin2,fill=-value))+facet_wrap(~name)+scale_fill_gradient2(na.value = "white")
    #
    #convert back value to the actual signal
    if (ref=="expected") {
      mat[,c("cv.group","phi"):=list(NULL,value*phihat.sd)]
      mat[,signal.old:=signal]
      mat[,signal:=exp(phi)]
      mat=mat[,.(name,bin1,bin2,phi,signal,signal.old,ncounts,value)]
      if(mat[,all(abs(signal-signal.old)<tol)]) break
    } else {
      mat[,c("cv.group","delta"):=list(NULL,value*deltahat.sd)]
      mat[,diffsig.old:=diffsig]
      mat[,diffsig:=exp(delta)]
      mat[,phi.ref:=(phihat.ref/sigmasq.ref + (phihat-delta)/sigmasq)/(1/sigmasq.ref+1/sigmasq)]
      mat=mat[,.(name,bin1,bin2,phi.ref,delta,diffsig,diffsig.old,ncounts,value)]
      if(mat[,all(abs(diffsig-diffsig.old)<tol)]) break
    }
  }
  mat[,ncounts:=NULL]
  mat = csnorm:::detect_binless_patches(mat)
  return(mat)
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
    p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=-(pmax(-2,pmin(2,value)))))+
      geom_raster(aes(begin2,begin1,fill=-(pmax(-2,pmin(2,value)))))+
      guides(fill=F)+scale_fill_gradient2()+facet_wrap(~name)+
      geom_polygon(aes(begin2,begin1,group=patchno),colour="black",fill=NA,data=a)
  }
  print(p)
}

#' Binless detection of significant interactions wrt expected
#' 
#' @param cs CSnorm object
#' @param resolution,group see
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to
#'   identify the input matrices.
#' @param ncores number of cores used for parallelization
#'   
#' @return the binned matrix with additional information relating to these 
#'   significant interactions
#' @export
#' 
#' @examples
detect_binless_interactions = function(cs, resolution, group, ncores=1, niter=10, tol=1e-3,
                                       cv.fold=10, cv.gridsize=30, verbose=T){
  #get CSmat object
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csm, type="binteractions", threshold=-1, ref="expected", raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat = csnorm_detect_binless(cs, resolution=resolution, group=group, ref="expected", ncores=ncores,
                              niter=niter, tol=tol, cv.fold=cv.fold, cv.gridsize=cv.gridsize, verbose=verbose)
  csi=new("CSinter", mat=mat, type="binteractions", threshold=-1, ref="expected")
  #store back
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
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
detect_binless_differences = function(cs, resolution, group, ref, niter=10, tol=1e-3, ncores=1,
                                      cv.fold=10, cv.gridsize=30, verbose=T){
  idx1=get_cs_binned_idx(cs, resolution, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, raise=T)
  csm=csb@grouped[[idx2]]
  if (get_cs_interaction_idx(csm, type="bdifferences", threshold=-1, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  mat = csnorm_detect_binless(cs, resolution=resolution, group=group, ref=ref, ncores=ncores,
                              niter=niter, tol=tol, cv.fold=cv.fold, cv.gridsize=cv.gridsize, verbose=verbose)
  csi=new("CSinter", mat=mat, type="bdifferences", threshold=-1, ref=ref)
  #add interaction to cs object
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}
  
