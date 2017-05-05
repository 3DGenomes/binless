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
  #diagonal
  chains=c(chains,list(c(sapply(chains,function(x){x[1]}),ntotal)))
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
                         alpha=0.2, inflate=2, tol.value=1e-6, maxsteps=100000) {
  #assume lambda1=0 and compute the fused lasso solution, centered around eCprime
  z=rep(0,tail(trails$breakpoints,n=1))
  u=rep(0,tail(trails$breakpoints,n=1))
  value = csnorm:::weighted_graphfl(valuehat, weight, trails$ntrails, trails$trails,
                                    trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.value/2, z, u) - eCprime
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
  submat[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, lambda1=0, lambda2=lambda2,
                                        eCprime=eCprime, tol.value=tol.value)]
  #get the number of patches and deduce degrees of freedom
  cl = csnorm:::build_patch_graph(submat, trails, tol.value=tol.value)$components
  submat[,patchno:=cl$membership]
  dof = submat[abs(value)+tol.value>lambda1,uniqueN(patchno)] #sparse fused lasso
  stopifnot(dof<=cl$no)
  #now soft-threshold the value around eCprime
  submat[,value:=sign(value)*pmax(abs(value)-lambda1, 0)]
  #compute BIC
  BIC = submat[,sum(weight*((valuehat-(value+eCprime))^2))+log(sum(ncounts))*dof]
  #compute mallow's Cp
  #Cp = submat[,sum(weight*((valuehat-(value+eCprime))^2 - 1))]+2*dof
  return(BIC)
}

#' cross-validate lambda1 and eCprime
#' 
#' @keywords internal
optimize_lambda1_eCprime = function(matg, trails, tol.val=1e-3, lambda2=0, positive=F,
                                    lambda1.min=0.05, constrained=T) {
  #compute values for lambda1=0 and eCprime=0
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, lambda2, 0, tol.value=tol.val)]
  matg[,value.ori:=value]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=factor(round(value,3))))+guides(fill=F))
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value<0.30)))
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, trails, tol.value=tol.val)$components
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1],size=.N),by=patchno][,.(N=.N,size=size[1]),keyby=value]
  if (patches[,.N]==1) #matrix is just one big bin
    return(data.table(eCprime=patches[,value],lambda1=lambda1.min,
                      BIC=matg[,sum(weight*((valuehat-(value+patches[,value]))^2))],dof=0))
  stopifnot(patches[,sum(N)]==cl$no)
  minval=patches[,min(value)]
  maxval=patches[,max(value)]
  valrange=maxval-minval
  #
  if (constrained==T) { #some patches must be zeroed to avoid degeneracy with diagonal decay fit
    forbidden.vals = matg[,.(max(value.ori)-min(value.ori)<=tol.val,value.ori[1]),by=diag.idx][V1==T,unique(V2)]
    if (positive==T) {
      lambda1.min = max(lambda1.min, (max(forbidden.vals)-minval)/2+tol.val)
    } else {
      lambda1.min = max(lambda1.min, (max(forbidden.vals)-min(forbidden.vals))/2+tol.val)
    }
  } else {
    forbidden.vals = c()
  }
  #
  obj = function(lambda1) {
    eCvals = c(patches[,value-lambda1],patches[,value+lambda1])
    if (positive==T) eCvals=eCvals[eCvals<=lambda1+minval]
    if (constrained==T) for (fv in forbidden.vals) eCvals = eCvals[abs(eCvals-fv)<=lambda1]
    if (length(eCvals)==0) 
      return(data.table(eCprime=0,lambda1=lambda1,BIC=.Machine$double.xmax,dof=NA))
    dt = foreach (eCprime=eCvals, .combine=rbind) %do% {
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
    dt[BIC==min(BIC)][lambda1==max(lambda1)][eCprime==min(eCprime)]
  }
  #dt=foreach(lambda1=seq(tol.val/2,valrange,length.out=50), .combine=rbind) %do% obj(lambda1)
  #dt=foreach(lambda1=seq(0,0.1,length.out=50), .combine=rbind) %do% obj(lambda1)
  #ggplot(dt)+geom_point(aes(lambda1,BIC,colour=ori))+geom_line(aes(lambda1,BIC,colour=ori))
  #ggplot(melt(dt,id.vars=c("lambda1","ori")))+geom_point(aes(lambda1,value,colour=ori))+
  #  geom_line(aes(lambda1,value,colour=ori))+facet_wrap(~variable,scales="free")
  #dt[,g:=dof*matg[,log(sum(ncounts))]]
  #dt[,h:=-2*matg[,.N]*log(lambda1)]
  #dt[,h:=-9*log(lambda1)+5*lambda1]
  #dt[,f:=BIC-g-h]
  #ggplot(dt)+geom_line(aes(lambda1,f,colour="f"))+geom_line(aes(lambda1,g,colour="g"))+geom_line(aes(lambda1,h,colour="h"))
  #ggplot(data.table(x=seq(0,valrange,l=1000))[,.(x,y=dgamma(x,rate=10,shape=5,log=T))])+geom_line(aes(x,y))
  #if (valrange <= 2*lambda1.min) return(as.list(obj(lambda1.min)))
  op=optimize(function(x){obj(10^(x))[,BIC]}, c(log10(max(lambda1.min,tol.val/2)),log10(valrange)), tol=tol.val)
  values=as.list(obj(10^(op$minimum)))
  #patches[,removed:=abs(value-values$eCprime)<=values$lambda1]
  return(values)
}

#' cross-validate lambda1 and eCprime, assuming positive=T and holding eCprime=lambda1+min(value)
#' 
#' @keywords internal
optimize_lambda1_eCprime_simplified = function(matg, trails, tol.val=1e-3, lambda2=0,
                                               lambda1.min=0.05, constrained=T) {
  #compute values for lambda1=0 and eCprime=0
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, lambda2, 0, tol.value=tol.val)]
  matg[,value.ori:=value]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=factor(round(value,3))))+guides(fill=F))
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value<0.30)))
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, trails, tol.value=tol.val)$components
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1],size=.N),by=patchno][,.(N=.N,size=size[1]),keyby=value]
  if (patches[,.N]==1) #matrix is just one big bin
    return(data.table(eCprime=patches[,value],lambda1=lambda1.min,
                      BIC=matg[,sum(weight*((valuehat-(value+patches[,value]))^2))],dof=0))
  stopifnot(patches[,sum(N)]==cl$no)
  minval=patches[,min(value)]
  maxval=patches[,max(value)]
  valrange=maxval-minval
  #
  if (constrained==T) { #some patches must be zeroed to avoid degeneracy with diagonal decay fit
    forbidden.vals = matg[,.(max(value.ori)-min(value.ori)<=tol.val,value.ori[1]),by=diag.idx][V1==T,unique(V2)]
    lambda1.min = max(lambda1.min, (max(forbidden.vals)-minval)/2 + tol.val)
  } else {
    forbidden.vals = c()
  }
  #
  obj = function(lambda1) {
    if (valrange<2*lambda1+tol.val) {
      eCprime=(maxval+minval)/2
    } else {
      eCprime=lambda1+minval-tol.val
    }
    if (constrained==T & any(abs(eCprime-forbidden.vals)>lambda1+tol.val))
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
  #dt=foreach(lambda1=seq(0,0.1,length.out=50), .combine=rbind) %do% obj(lambda1)
  #ggplot(dt)+geom_point(aes(lambda1,BIC,colour=ori))+geom_line(aes(lambda1,BIC,colour=ori))
  #ggplot(melt(dt,id.vars=c("lambda1","ori")))+geom_point(aes(lambda1,value,colour=ori))+
  #  geom_line(aes(lambda1,value,colour=ori))+facet_wrap(~variable,scales="free")
  #dt[,g:=dof*matg[,log(sum(ncounts))]]
  #dt[,h:=-2*matg[,.N]*log(lambda1)]
  #dt[,h:=-9*log(lambda1)+5*lambda1]
  #dt[,f:=BIC-g-h]
  #ggplot(dt)+geom_line(aes(lambda1,f,colour="f"))+geom_line(aes(lambda1,g,colour="g"))+geom_line(aes(lambda1,h,colour="h"))
  #ggplot(data.table(x=seq(0,valrange,l=1000))[,.(x,y=dgamma(x,rate=10,shape=5,log=T))])+geom_line(aes(x,y))
  if (valrange <= 2*lambda1.min) return(as.list(obj(lambda1.min)))
  op=optimize(function(x){obj(10^(x))[,BIC]}, c(log10(max(lambda1.min,tol.val/2)),log10(valrange)), tol=tol.val)
  values=as.list(obj(10^(op$minimum)))
  #patches[,removed:=abs(value-values$eCprime)<=values$lambda1]
  return(values)
}

#' cross-validate lambda1 and assume eCprime=0
#' 
#' @keywords internal
optimize_lambda1_only = function(matg, trails, tol.val=1e-3, lambda2=0, positive=F,
                                 lambda1.min=0.05, constrained=T) {
  #compute values for lambda1=0 and eCprime=0
  matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, lambda2, 0, tol.value=tol.val)]
  matg[,value.ori:=value]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=factor(round(value,3))))+guides(fill=F))
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value<0.30)))
  #get the number of patches to compute degrees of freedom
  cl = csnorm:::build_patch_graph(matg, trails, tol.value=tol.val)$components
  matg[,patchno:=cl$membership]
  patches = matg[,.(value=value[1],size=.N),by=patchno][,.(N=.N,size=size[1]),keyby=value]
  if (patches[,.N]==1) if (abs(patches[,value])<=lambda1.min) #matrix is just one big bin
    return(data.table(lambda1=lambda1.min,
                      BIC=matg[,sum(weight*(valuehat-value)^2)],dof=0))
  stopifnot(patches[,sum(N)]==cl$no)
  minval=patches[,min(abs(value))]
  maxval=patches[,max(abs(value))]
  if (positive==T) lambda1.min=max(lambda1.min, minval)
  #
  if (constrained==T) { #some patches must be zeroed to avoid degeneracy with diagonal decay fit
    forbidden.vals = matg[,.(max(value.ori)-min(value.ori)<=tol.val,value.ori[1]),by=diag.idx][V1==T,unique(V2)]
    lambda1.min = max(lambda1.min, forbidden.vals)
  } else {
    forbidden.vals = c()
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
    data.table(lambda1=lambda1,BIC=BIC,dof=dof)
  }
  #dt=foreach(lambda1=seq(lambda1.min,maxval,length.out=50), .combine=rbind) %do% obj(lambda1)
  #ggplot(dt)+geom_point(aes(lambda1,BIC))+geom_line(aes(lambda1,BIC))
  #dt[,g:=dof*matg[,log(sum(ncounts))]]
  #dt[,h:=-2*matg[,.N]*log(lambda1)]
  #dt[,h:=-9*log(lambda1)+5*lambda1]
  #dt[,f:=BIC-g-h]
  #ggplot(dt)+geom_line(aes(lambda1,f,colour="f"))+geom_line(aes(lambda1,g,colour="g"))+geom_line(aes(lambda1,h,colour="h"))
  #ggplot(data.table(x=seq(0,valrange,l=1000))[,.(x,y=dgamma(x,rate=10,shape=5,log=T))])+geom_line(aes(x,y))
  op=optimize(function(x){obj(10^(x))[,BIC]}, c(log10(max(lambda1.min,tol.val/2)),log10(maxval)), tol=tol.val)
  return(as.list(obj(10^(op$minimum))))
}

#' cross-validate lambda2
#' @keywords internal
optimize_lambda2 = function(matg, trails, tol.val=1e-3, lambda2.max=1000) {
  obj = function(x){csnorm:::gfl_BIC(matg, trails, lambda1=0, lambda2=10^(x), eCprime=0, tol.value = tol.val)}
  minlambda=tol.val*10
  maxlambda=lambda2.max
  #save(matg,trails,tol,lambda1,eCprime,file="debug_lambda2.RData")
  #cat("*** maxlambda ",maxlambda," (range=",matg[,max(value)-min(value)],")\n")
  #shrink maximum lambda, in case initial guess is too big
  repeat {
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, maxlambda, 0, tol.value=tol.val)]
    #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+geom_raster(aes(bin2,bin1,fill=valuehat)))
    if (matg[,max(value)-min(value)] > tol.val) {
      maxlambda = maxlambda*2
      break
    }
    maxlambda = maxlambda/2
    #dof = csnorm:::get_gfl_degrees_of_freedom(matg, trails, tol.val=tol.val)
    #cat("shrink maxlambda ",maxlambda," (range=",matg[,max(value)-min(value)],", dof=",dof,")\n")
  }
  #cat("optimize lambda2 between ",minlambda," and ",maxlambda,"\n")
  #dt = foreach (lam=seq(log10(minlambda),log10(maxlambda),l=100),.combine=rbind) %dopar% data.table(x=lam,y=obj(lam))
  #ggplot(dt)+geom_point(aes(10^(x),y))#+scale_x_log10()+scale_y_log10()
  op=optimize(obj, c(log10(minlambda),log10(maxlambda)), tol=tol.val)
  lambda2=10^op$minimum
  if (lambda2==maxlambda) cat("   Warning: lambda2 hit upper boundary.\n")
  if (lambda2 <= tol.val*10) {
    cat("   Warning: lambda2 too close to lower boundary.")
    lambda2=0
  }
  return(lambda2)
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
#' @param ncores integer (default 1).
#'   
#'   finds optimal lambda1, lambda2 and eC using BIC.
#' @keywords internal
csnorm_fused_lasso = function(matg, trails, positive, fixed, constrained, simplified, tol.val=1e-3,
                              verbose=T, ncores=1) {
  groupname=matg[,name[1]]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=valuehat))+scale_fill_gradient2())
  lambda2 = csnorm:::optimize_lambda2(matg, trails, tol.val = tol.val)
  #get best lambda1 and set eCprime to lower bound
  if (fixed==F) {
    if (simplified==T) {
      if (positive!=T) stop("Cannot have positive!=T and simplified==T")
      vals = csnorm:::optimize_lambda1_eCprime_simplified(matg, trails, tol.val=tol.val, lambda2=lambda2,
                                                          constrained=constrained)
    } else {
      vals = csnorm:::optimize_lambda1_eCprime(matg, trails, tol.val=tol.val, lambda2=lambda2,
                                               constrained=constrained, positive=positive)
    }
    eCprime=vals$eCprime
  } else {
    vals = csnorm:::optimize_lambda1_only(matg, trails, tol.val=tol.val, lambda2=lambda2, positive=positive,
                                          constrained=constrained)
    vals$eCprime=0
  }
  vals$lambda2=lambda2
  vals$name=groupname
  #matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, 0, vals$lambda2, 0, tol.value = tol.val)]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=abs(value-vals$eCprime)<=vals$lambda1)))
  #matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, vals$lambda1, vals$lambda2, vals$eCprime, tol.value = tol.val)]
  #print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())
  return(as.data.table(vals))
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
  return(list(mat=mat,trails=trails,eCprime=names[,.(name=groupname,eCprime=0)]))
}

#' Build grouped difference matrix using normalization data if available
#' @keywords internal
prepare_difference_matrix = function(cs, names, resolution, ref) {
  stuff = csnorm:::prepare_signal_matrix(cs, names, resolution)
  mat = foreach(n=names[groupname!=ref,unique(groupname)],.combine=rbind) %do%
    merge(stuff$mat[name==n],stuff$mat[name==ref,.(bin1,bin2,phi1=phi)],all=T,by=c("bin1","bin2"))
  mat[,c("phi.ref","delta"):=list((phi+phi1)/2,(phi-phi1)/2)]
  mat[,c("phi","phi1"):=NULL]
  return(list(mat=mat,trails=stuff$trails,eCprime=names[groupname!=ref,.(name=groupname,eCprime=0)]))
}


#' compute input to fused lasso
#' @keywords internal
csnorm_compute_raw_signal = function(cts, dispersion, mat, eCprime) {
  mat=mat[,.(name,bin1,bin2,phi)]
  cts.cp = mat[eCprime[cts[,.(name,bin1,bin2,count,lmu.nosig,weight,var)],,on="name"],,on=c("name","bin1","bin2")]
  cts.cp[,mu:=exp(lmu.nosig+phi+eCprime)]
  cts.cp[,c("z","var"):=list(count/mu-1, (1/mu+1/dispersion))]
  mat = cts.cp[,.(phihat=weighted.mean(z+phi, weight/var),
                  phihat.var=2/sum(weight/var),
                  ncounts=sum(weight)),keyby=c("name","bin1","bin2")][mat,,on=c("name","bin1","bin2")]
  mat[is.na(phihat),c("phihat","phihat.var","ncounts"):=list(1,Inf,0)] #bins with no detectable counts
  mat[,c("valuehat","weight"):=list(phihat,1/phihat.var)]
  setkey(mat,name,bin1,bin2)
  return(mat)
}

#' compute input to fused lasso
#' @keywords internal
csnorm_compute_raw_differential = function(cts, dispersion, mat, eCprime, ref) {
  ctsref = foreach(n=cts[name!=ref,unique(name)],.combine=rbind) %do%
    cts[name==ref,.(name=n,bin1,bin2,count,lmu.nosig,weight,var)]
  ctsoth=cts[name!=ref,.(name,bin1,bin2,count,lmu.nosig,weight,var)]
  mat=mat[,.(name,bin1,bin2,phi.ref,delta)]
  #
  mat.ref = csnorm:::csnorm_compute_raw_signal(ctsref,dispersion,mat[,.(name,bin1,bin2,phi=phi.ref)], eCprime[,.(name,eCprime=0)])
  mat.ref[,c("valuehat","weight"):=NULL]
  setnames(mat.ref,c("phihat","phihat.var","ncounts","phi"),c("phihat.ref","phihat.var.ref","ncounts.ref","phi.ref"))
  mat.oth = csnorm:::csnorm_compute_raw_signal(ctsoth,dispersion,mat[,.(name,bin1,bin2,phi=phi.ref+delta)],eCprime)
  stopifnot(mat.oth[,.N]==mat.ref[,.N])
  #
  mat=merge(mat.oth,mat.ref)
  mat[,c("deltahat","deltahat.var","ncounts"):=list(phihat-phihat.ref,phihat.var+phihat.var.ref,ncounts+ncounts.ref)]
  mat[,c("valuehat","weight","ncounts.ref","delta"):=list(deltahat,1/deltahat.var,NULL,phi-phi.ref)]
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
detect_binless_interactions = function(cs, resolution, group, ncores=1, niter=10, tol.val=1e-5, verbose=T){
  if (verbose==T) cat("Binless interaction detection with resolution=",resolution," and group=",group,"\n")
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csg, type="binteractions", threshold=-1, ref="expected", raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  #
  ### prepare signal estimation
  if (verbose==T) cat("  Prepare for signal estimation\n")
  stuff = csnorm:::prepare_signal_matrix(cs, csg@names, resolution)
  mat=stuff$mat
  trails=stuff$trails
  eCprime=stuff$eCprime
  ### main loop
  registerDoParallel(cores=ncores)
  for (step in 1:niter) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_signal(csg@cts, csg@par$alpha, mat, eCprime)
    #ggplot(mat)+geom_raster(aes(bin2,bin1,fill=phi))+geom_raster(aes(bin1,bin2,fill=phihat))+scale_fill_gradient2()+facet_wrap(~name)
    #
    #perform fused lasso on signal, at fixed offset
    if (verbose==T) cat("  Fused lasso\n")
    groupnames=mat[,unique(name)]
    params = foreach(g=groupnames, .combine=rbind) %dopar%
      csnorm:::csnorm_fused_lasso(mat[name==g], trails, fixed=F, positive=T, constrained=F, simplified=F, tol.val=tol.val,
                                  ncores=ncores, verbose=verbose)
    #display param info
    if (verbose==T)
      for (i in 1:params[,.N])
        cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
            " eCprime=",params[i,eCprime],"\n")
    #compute matrix at new params
    #save(mat,params,file=paste0("mat_step_",step,".RData"))
    mat = foreach (g=groupnames, .combine=rbind) %dopar% {
      p=params[name==g]
      matg=mat[name==g]
      matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, p$lambda1, p$lambda2, p$eCprime, tol.value = tol.val)]
      matg
    }
    eCprime = eCprime[params[,.(name,eCprime2=eCprime)],,on="name"][,.(name,eCprime=eCprime+eCprime2)]
    #convert back value to the actual signal
    mat[,phi.old:=phi]
    mat[,phi:=value]
    #
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("sig_step_",step,"_value.png"), width=10, height=8)
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("sig_step_",step,"_weight.png"), width=10, height=8)
    #
    #check convergence
    if(mat[,all(abs(phi-phi.old)<tol.val)] & eCprime[,all(abs(eCprime)<tol.val)]) break
  }
  #
  if (verbose==T) cat(" Detect patches\n")
  mat = csnorm:::detect_binless_patches(mat, trails, tol.value=tol.val)
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
detect_binless_differences = function(cs, resolution, group, ref, ncores=1, niter=10, tol.val=1e-3, verbose=T){
  ### get CSgroup object
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  if (get_cs_interaction_idx(csg, type="bdifferences", threshold=-1, ref=ref, raise=F)>0)
    stop("Refusing to overwrite this already detected interaction")
  if (is.character(ref)) ref=csg@names[as.character(groupname)==ref,unique(groupname)]
  if (verbose==T) cat("Binless difference detection with resolution=",resolution,
                        " group=",group," and ref=",as.character(ref),"\n")
  if (verbose==T) cat("  Prepare for difference estimation\n")
  stuff = csnorm:::prepare_difference_matrix(cs, csg@names, resolution, ref)
  mat=stuff$mat
  trails=stuff$trails
  eCprime=stuff$eCprime
  ### main loop
  registerDoParallel(cores=ncores)
  for (step in 1:niter) {
    if (verbose==T) cat(" Main loop, step ",step,"\n")
    if (verbose==T) cat("  Estimate raw signal\n")
    mat = csnorm:::csnorm_compute_raw_differential(csg@cts, csg@par$alpha, mat, eCprime, ref)
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=deltahat))+facet_wrap(~name)+scale_fill_gradient2()
    #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi.ref))+geom_raster(aes(bin2,bin1,fill=phihat))+facet_wrap(~name)+scale_fill_gradient2()
    #
    #perform fused lasso on signal
    if (verbose==T) cat("  Fused lasso\n")
    groupnames=mat[,unique(name)]
    params = foreach(g=groupnames, .combine=rbind) %dopar%
      csnorm:::csnorm_fused_lasso(mat[name==g], trails, fixed=F, positive=F, constrained=F, simplified=F,
                                  tol.val=tol.val, ncores=ncores, verbose=verbose)
    #display param info
    if (verbose==T)
      for (i in 1:params[,.N])
        cat("  ",params[i,name]," : lambda1=",params[i,lambda1]," lambda2=",params[i,lambda2],
            " eC'=",params[i,eCprime],"\n")
    #compute matrix at new params
    #save(mat,params,file=paste0("dmat_step_",step,".RData"))
    mat = foreach (g=groupnames, .combine=rbind) %dopar% {
      p=params[name==g]
      matg=mat[name==g]
      matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, p$lambda1, p$lambda2, p$eCprime, tol.value = tol.val)]
      matg
    }
    eCprime = eCprime[params[,.(name,eCprime2=eCprime)],,on="name"][,.(name,eCprime=eCprime+eCprime2)]
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("diff_step_",step,"_value.png"), width=10, height=8)
    #p=ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))+scale_fill_gradient2()+facet_wrap(~name)
    #ggsave(p,filename = paste0("diff_step_",step,"_weight.png"), width=10, height=8)
    #
    #convert back value to the actual signal
    mat[,delta.old:=delta]
    mat[,delta:=value]
    mat[,phi.ref:=(phihat.ref/phihat.var.ref + (phihat-delta)/phihat.var)/(1/phihat.var.ref+1/phihat.var)]
    if(mat[,all(abs(delta-delta.old)<tol.val)]  & eCprime[,all(abs(eCprime)<tol.val)]) break
  }
  if (verbose==T) cat(" Detect patches\n")
  mat = csnorm:::detect_binless_patches(mat, trails, tol.value=tol.val)
  csi=new("CSinter", mat=mat, type="bdifferences", threshold=-1, ref=as.character(ref))
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
