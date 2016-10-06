#' @include csnorm.R
NULL

#' Convert a sparse counts data.table to a dense one by adding rows with zero counts
#'
#' @param counts,biases data.tables as returned by \code{\link{prepare_for_sparse_cs_norm}}
#' @param biases2 data.table of biases for id2 column of counts. If NULL (default), use that of biases
#'
#' @return a counts data.table with zeros filled according to cut sites provided in biases (and biases2 if available)
#' @keywords internal
#' @export
#' @section Warning:
#' Memory-intensive
#' @examples
fill_zeros = function(counts,biases,biases2=NULL,circularize=-1L) {
  if (is.null(biases2)) biases2=biases
  runname=biases[,unique(name)]
  if (length(runname)>1) {
    foreach (i=runname, .combine=rbind) %do%
      fill_zeros(counts[name==i],biases[name==i],biases2[name==i],circularize=circularize)
  } else {
    if (biases[,.N]==0 | biases2[,.N]==0) return(data.table())
    newcounts=CJ(biases[,paste(id,pos)],biases2[,paste(id,pos)])
    newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
    newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
    newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=
                list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
    newcounts=newcounts[pos1<pos2]
    setkey(newcounts, id1, id2, pos1, pos2)
    setkey(counts, id1, id2, pos1, pos2)
    newcounts=counts[newcounts]
    newcounts[is.na(contact.close),contact.close:=0]
    newcounts[is.na(contact.far),contact.far:=0]
    newcounts[is.na(contact.up),contact.up:=0]
    newcounts[is.na(contact.down),contact.down:=0]
    newcounts[is.na(name),name:=runname]
    if (circularize>0) {
      newcounts[,distance:=pmin(abs(pos2-pos1), circularize+1-abs(pos2-pos1))]
    } else {
      newcounts[,distance:=abs(pos2-pos1)]
    }
    newcounts
  }
}

#' Get subset of counts and biases data.table
#'
#' @param counts,biases as returned by \code{\link{prepare_for_sparse_cs_norm}}
#' @param begin1,end1 integers. Genomic location of begin/end of subset to take 
#' @param begin2,end2 integers. If provided, take extradiagonal section 
#' @param fill.zeros whether to return a dense matrix (default TRUE)
#' @param integer. Set this to the size of the chromosome if it is
#'   circular, otherwise leave as-is (default is -1)
#'   
#' @return counts and biases data.tables, in a list
#' @export
#' @section Warning:
#' Memory-intensive
#' @examples
get_cs_subset = function(counts, biases, begin1, end1, begin2=NULL, end2=NULL, fill.zeros=T, circularize=-1) {
  stopifnot(end1>=begin1)
  if (is.null(begin2)) begin2=begin1
  if (is.null(end2)) end2=end1
  stopifnot(begin2>=begin1)
  #
  biases.1=biases[pos>=begin1&pos<end1]
  beginrange1=biases.1[1,id]
  endrange1=biases.1[.N,id]
  biases.1[,id:=id-beginrange1+1]
  #
  biases.2=biases[pos>=begin2&pos<end2]
  beginrange2=biases.2[1,id]
  endrange2=biases.2[.N,id]
  biases.2[,id:=id-beginrange2+1]
  #
  counts.local=counts[pos1>=begin1&pos1<end1&pos2>=begin2&pos2<end2]
  counts.local[,c("id1","id2"):=list(id1-beginrange1+1,id2-beginrange2+1)]
  #fill zeros
  if (fill.zeros==T) counts.local=fill_zeros(counts=counts.local,biases=biases.1,biases2=biases.2, circularize=circularize)
  return(list(counts=counts.local,biases1=biases.1,biases2=biases.2,
              beginrange1=beginrange1, endrange1=endrange1, beginrange2=beginrange2, endrange2=endrange2))
}

#' Single-cpu simplified initial guess
#' @keywords internal
#' 
csnorm_simplified_guess = function(biases, counts, design, lambda, dmin, dmax, bf_per_kb=1, bf_per_decade=5, groups=10,
                                   iter=10000, dispersion=10) {
  for (n in biases[,unique(name)]) {
    stopifnot(groups<=biases[name==n,.N-1])
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #collect all counts on left/right side and put into quantile groups
  cts=rbind(counts[,.(id=id1,ldist=log(distance),R=(contact.close+contact.down),L=(contact.far+contact.up),others=2)],
           counts[,.(id=id2,ldist=log(distance),R=(contact.far+contact.down),L=(contact.close+contact.up),others=2)])
  setkey(cts,id)
  cts[,cbin:=ntile(L+R,groups),by=id]
  ctsl=dcast(cts[,.(id,cbin,L)], id~cbin, value.var="L", fun.aggregate=sum)
  ctsl[,id:=NULL]
  stopifnot(dim(ctsl)==c(biases[,.N],groups))
  ctsr=dcast(cts[,.(id,cbin,R)], id~cbin, value.var="R", fun.aggregate=sum)
  ctsr[,id:=NULL]
  stopifnot(dim(ctsr)==c(biases[,.N],groups))
  ctso=dcast(cts[,.(id,cbin,others)], id~cbin, value.var="others", fun.aggregate=sum)
  ctso[,id:=NULL]
  stopifnot(dim(ctso)==c(biases[,.N],groups))
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)], XB=as.array(design[,genomic]),
            Krow=Krow, SD=biases[,.N], bbegin=bbegin,
            cutsitesD=biases[,pos], rejoined=biases[,rejoined],
            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
            G=groups, counts_sum_left=ctsl, counts_sum_right=ctsr, log_decay_sum=log(ctso), lambda=lambda, alpha=dispersion)
  op=optimizing(csnorm:::stanmodels$simplified_guess, data=data, as_vector=F, hessian=F, iter=iter, verbose=T, init=0, init_alpha=1e-9)
  #add diagonal decay inits
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  op$par=c(list(lambda_nu=array(lambda,dim=data$Biases), lambda_delta=array(lambda,dim=data$Biases),
                beta_diag=beta_diag, lambda_diag=array(1,dim=Decays), log_decay=rep(0,counts[,.N])), alpha=dispersion,
           op$par[c("eC","eRJ","eDE","beta_nu","beta_delta","log_nu","log_delta")])
  return(op)
}

#' Single-cpu simplified fitting for nu and delta
#' @keywords internal
#' 
csnorm_simplified_decay = function(biases, counts, design, log_nu, log_delta, dmin, dmax, 
                                   bf_per_decade=5, bins_per_bf=10, groups=10,
                                   iter=10000, verbose=T, init=0, dispersion=10, ...) {
  for (n in biases[,unique(name)]) {
    stopifnot(groups<=biases[name==n,.N-1])
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #add bias informations to counts
  setkey(counts, id1, id2)
  csub=copy(counts)
  bsub=biases[,.(id)]
  bsub[,c("nu","delta"):=list(exp(log_nu),exp(log_delta))]
  csub=merge(bsub[,.(id1=id,nu,delta)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,nu,delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  #bin distances
  stepsz=1/(bins_per_bf*bf_per_decade)
  dbins=10**seq(log10(dmin),log10(dmax)+stepsz,stepsz)
  csub[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #collect all counts in these bins and put into quantile groups
  csub[,c("ldist","count","others"):=list(
    log(distance), contact.far+contact.down+contact.close+contact.up, nu1*nu2*(delta1+1/delta1)*(delta2+1/delta2))]
  csub[,cbin:=ntile(count,groups),by=dbin]
  csd = csub[,.(mdist=exp(mean(ldist)), count=sum(count), others=sum(others), weight=4*.N), keyby=c("name", "dbin","cbin")]
  #run optimization
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  data=list(Dsets=design[,.N], Decays=design[,uniqueN(decay)], XD=as.array(design[,decay]),
            Kdiag=Kdiag, dmin=dmin, dmax=dmax, N=csd[,.N], cbegin=cbegin,
            counts_sum=csd[,count], weight=csd[,weight], dist=csd[,mdist], log_genomic_sum=csd[,log(others)],
            alpha=dispersion)
  op=optimizing(stanmodels$simplified_decay, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init,
                init_alpha=1e-9, tol_rel_grad=0, tol_rel_obj=1e3, ...)
  #make nice decay data table
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  #rewrite log_decay as if it was calculated for each count
  csd[,log_decay:=op$par$log_decay]
  a=csd[csub,.(id1,id2,pos1,pos2,distance,log_decay),on=key(csd)]
  setkeyv(a,key(counts))
  op$par$log_decay=a[,log_decay]
  return(op)
}

#' Single-cpu simplified fitting for nu and delta
#' @keywords internal
#' 
csnorm_simplified_genomic = function(biases, counts, design, log_decay, log_nu, log_delta, bf_per_kb=1, groups=10,
                                     iter=10000, verbose=T, init=0, dispersion=10, ...) {
  for (n in biases[,unique(name)]) {
    stopifnot(groups<=biases[name==n,.N-1])
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #add bias informations to counts
  csub=copy(counts)
  csub[,decay:=exp(log_decay)]
  bsub=biases[,.(id)]
  bsub[,c("nu","delta"):=list(exp(log_nu),exp(log_delta))]
  csub=merge(bsub[,.(id1=id,nu,delta)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,nu,delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  #collect all counts on left/right side and put into quantile groups
  cts=rbind(csub[,.(id=id1,ldist=log(distance),R=(contact.close+contact.down),L=(contact.far+contact.up),others=decay*nu2*(delta2+1/delta2))],
            csub[,.(id=id2,ldist=log(distance),R=(contact.far+contact.down),L=(contact.close+contact.up),others=decay*nu1*(delta1+1/delta1))])
  setkey(cts,id)
  cts[,cbin:=ntile(L+R,groups),by=id]
  ctsl=dcast(cts[,.(id,cbin,L)], id~cbin, value.var="L", fun.aggregate=sum)
  ctsl[,id:=NULL]
  stopifnot(dim(ctsl)==c(biases[,.N],groups))
  ctsr=dcast(cts[,.(id,cbin,R)], id~cbin, value.var="R", fun.aggregate=sum)
  ctsr[,id:=NULL]
  stopifnot(dim(ctsr)==c(biases[,.N],groups))
  ctso=dcast(cts[,.(id,cbin,others)], id~cbin, value.var="others", fun.aggregate=sum)
  ctso[,id:=NULL]
  stopifnot(dim(ctso)==c(biases[,.N],groups))
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)], XB=as.array(design[,genomic]),
            Krow=Krow, SD=biases[,.N], bbegin=bbegin,
            cutsitesD=biases[,pos], rejoined=biases[,rejoined],
            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
            G=groups, counts_sum_left=ctsl, counts_sum_right=ctsr, log_decay_sum=log(ctso),
            alpha=dispersion)
  optimizing(stanmodels$simplified_genomic, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose,
             init=init, init_alpha=1e-9, ...)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_simplified_dispersion = function(biases, counts, design, dmin, dmax, beta_nu, beta_delta, beta_diag,
                                        lambda_nu, lambda_delta, lambda_diag,
                                        bf_per_kb = 1, bf_per_decade=5, iter = 10000,
                                        verbose=T, init=0, weight=array(1,dim=design[,.N]),...) {
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=design[,.N], Biases=design[,uniqueN(genomic)], Decays=design[,uniqueN(decay)],
               XB=as.array(design[,genomic]), XD=as.array(design[,decay]),
               Krow=Krow, SD=biases[,.N], bbegin=bbegin,
               cutsitesD=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cbegin=cbegin,
               cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight), beta_nu=beta_nu, beta_delta=beta_delta, beta_diag=beta_diag,
               lambda_nu=lambda_nu, lambda_delta=lambda_delta, lambda_diag=lambda_diag)
  op=optimizing(stanmodels$simplified_dispersion, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose,
                init=init, init_alpha=1e-9, ...)
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  return(op)
}

#' Single-cpu simplified initial guess
#' @keywords internal
#' 
csnorm_gauss_guess = function(biases, counts, design, lambda, dmin, dmax, bf_per_kb=1, bf_per_decade=5,
                                   iter=10000, dispersion=10) {
  for (n in biases[,unique(name)]) {
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #collect all counts on left/right side. No need to put name since IDs are unique and increasing
  cts=rbind(counts[,.(id=id1,R=(contact.close+contact.down),L=(contact.far+contact.up),others=2)],
            counts[,.(id=id2,R=(contact.far+contact.down),L=(contact.close+contact.up),others=2)])
  cts=cts[,.(ctsl=sum(L),ctsr=sum(R),ctso=sum(others)),keyby=id]
  stopifnot(cts[,.N]==biases[,.N])
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)], XB=as.array(design[,genomic]),
            Krow=Krow, SD=biases[,.N], bbegin=bbegin,
            cutsitesD=biases[,pos], rejoined=biases[,rejoined],
            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
            counts_sum_left=cts[,ctsl], counts_sum_right=cts[,ctsr], log_decay_sum=cts[,log(ctso)],
            lambda=lambda, alpha=dispersion)
  op=optimizing(csnorm:::stanmodels$gauss_guess, data=data, as_vector=F, hessian=F, iter=iter,
                verbose=T, init=0, init_alpha=1e-9)
  #add diagonal decay inits
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  op$par=c(list(lambda_nu=array(lambda,dim=data$Biases), lambda_delta=array(lambda,dim=data$Biases),
                beta_diag=beta_diag, lambda_diag=array(1,dim=Decays), log_decay=rep(0,counts[,.N]),
                alpha=dispersion),
           op$par[c("eC","eRJ","eDE","beta_nu","beta_delta","log_nu","log_delta")])
  return(op)
}

#' Single-cpu simplified fitting for nu and delta
#' @keywords internal
#' 
csnorm_gauss_decay = function(biases, counts, design, log_nu, log_delta, dmin, dmax, dispersion, lambda_diag,
                              bf_per_decade=5, bins_per_bf=10, iter=10000, verbose=T, init=0, ...) {
  for (n in biases[,unique(name)]) {
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #add bias informations to counts
  setkey(counts, id1, id2)
  csub=copy(counts)
  bsub=biases[,.(id)]
  bsub[,c("nu","delta"):=list(exp(log_nu),exp(log_delta))]
  csub=merge(bsub[,.(id1=id,nu,delta)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,nu,delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  #bin distances
  stepsz=1/(bins_per_bf*bf_per_decade)
  dbins=10**seq(log10(dmin),log10(dmax)+stepsz,stepsz)
  csub[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #collect all counts in these bins
  csub[,c("ldist","count","others"):=list(
    log(distance), contact.far+contact.down+contact.close+contact.up, nu1*nu2*(delta1+1/delta1)*(delta2+1/delta2))]
  csd = csub[,.(mdist=exp(mean(ldist)), count=sum(count), others=sum(others), weight=4*.N),
             keyby=c("name", "dbin")]
  #run optimization
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  data=list(Dsets=design[,.N], Decays=design[,uniqueN(decay)], XD=as.array(design[,decay]),
            Kdiag=Kdiag, dmin=dmin, dmax=dmax, N=csd[,.N], cbegin=cbegin,
            counts_sum=csd[,count], weight=csd[,weight], dist=csd[,mdist], log_genomic_sum=csd[,log(others)],
            alpha=dispersion, lambda_diag=lambda_diag)
  op=optimizing(stanmodels$gauss_decay, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init,
                init_alpha=1e-9, tol_rel_grad=0, tol_rel_obj=1e3, ...)
  #make nice decay data table
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  #rewrite log_decay as if it was calculated for each count
  csd[,log_decay:=op$par$log_decay]
  a=csd[csub,.(id1,id2,pos1,pos2,distance,log_decay),on=key(csd)]
  setkeyv(a,key(counts))
  op$par$log_decay=a[,log_decay]
  return(op)
}

#' Single-cpu simplified fitting for nu and delta
#' @keywords internal
#' 
csnorm_gauss_genomic = function(biases, counts, design, log_decay, log_nu, log_delta, dispersion,
                                lambda_nu, lambda_delta, bf_per_kb=1, iter=10000, verbose=T, init=0, ...) {
  for (n in biases[,unique(name)]) {
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #add bias informations to counts
  csub=copy(counts)
  csub[,decay:=exp(log_decay)]
  bsub=biases[,.(id)]
  bsub[,c("nu","delta"):=list(exp(log_nu),exp(log_delta))]
  csub=merge(bsub[,.(id1=id,nu,delta)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,nu,delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  #collect all counts on left/right side
  cts=rbind(csub[,.(id=id1,ldist=log(distance),R=(contact.close+contact.down),L=(contact.far+contact.up),others=decay*nu2*(delta2+1/delta2))],
            csub[,.(id=id2,ldist=log(distance),R=(contact.far+contact.down),L=(contact.close+contact.up),others=decay*nu1*(delta1+1/delta1))])
  cts=cts[,.(ctsl=sum(L),ctsr=sum(R),ctso=sum(others)),keyby=id]
  stopifnot(cts[,.N]==biases[,.N])
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)], XB=as.array(design[,genomic]),
            Krow=Krow, SD=biases[,.N], bbegin=bbegin,
            cutsitesD=biases[,pos], rejoined=biases[,rejoined],
            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
            counts_sum_left=cts[,ctsl], counts_sum_right=cts[,ctsr], log_decay_sum=cts[,log(ctso)],
            alpha=dispersion, lambda_nu=lambda_nu, lambda_delta=lambda_delta)
  optimizing(stanmodels$gauss_genomic, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose,
             init=init, init_alpha=1e-9, ...)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_gauss_dispersion = function(biases, counts, design, dmin, dmax, beta_nu, beta_delta, beta_diag,
                                        bf_per_kb = 1, bf_per_decade=5, iter = 10000,
                                        verbose=T, init=0, weight=array(1,dim=design[,.N]),...) {
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=design[,.N], Biases=design[,uniqueN(genomic)], Decays=design[,uniqueN(decay)],
               XB=as.array(design[,genomic]), XD=as.array(design[,decay]),
               Krow=Krow, SD=biases[,.N], bbegin=bbegin,
               cutsitesD=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cbegin=cbegin,
               cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight), beta_nu=beta_nu, beta_delta=beta_delta, beta_diag=beta_diag)
  op=optimizing(stanmodels$gauss_dispersion, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose,
                init=init, init_alpha=1e-9, ...)
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  return(op)
}


#' Single-cpu fitting
#' @keywords internal
#' 
csnorm_fit = function(biases, counts, design, dmin, dmax, bf_per_kb=1, bf_per_decade=5, iter=10000,
                      verbose=T, init=0, weight=array(1,dim=design[,.N]), ...) {
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=design[,.N], Biases=design[,uniqueN(genomic)], Decays=design[,uniqueN(decay)],
               XB=as.array(design[,genomic]), XD=as.array(design[,decay]),
               Krow=Krow, SD=biases[,.N], bbegin=bbegin,
               cutsitesD=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cbegin=cbegin,
               cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight))
  op=optimizing(stanmodels$fit, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  return(op)
}

#' Single-cpu fitting, only diagonal decay
#' @keywords internal
#' 
csnorm_fit_extradiag = function(biases1, biases2, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, weight=1, init=0, ...) {
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( S1=biases1[,.N], S2=biases2[,.N], Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far], counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               log_nu1=biases1[,log_nu], log_nu2=biases2[,log_nu], log_delta1=biases1[,log_delta], log_delta2=biases2[,log_delta],
               weight=weight)
  optimizing(stanmodels$fit_extradiag, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

#' Predict expected values for each count given optimized model parameters
#' 
#' @param cs CSnorm object
#' @param counts
#' @param verbose
#'   
#' @section Warning: Do not call this function to compute a binned matrix. It
#'   should (almost) never be used and might be removed in the future.
#'   
#' @return a stan optimization output with the predictions.
#' @keywords internal
#' @export
#' 
#' @examples
csnorm_predict_all = function(cs, counts, verbose=T) {
  biases=cs@biases
  par=cs@par
  dmin=cs@settings$dmin
  dmax=cs@settings$dmax
  Kdiag=round((log10(dmax)-log10(dmin))*cs@settings$bf_per_decade)
  ncounts=counts[,.N]
  if (counts[,.N]==1) {
    cbegin=c(1,2)
  } else {
    cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  }
  design=cs@design
  data = list( Dsets=design[,.N], Decays=design[,uniqueN(decay)], XD=as.array(design[,decay]),
               Kdiag=Kdiag, SD=biases[,.N], cutsitesD=biases[,pos], dmin=dmin, dmax=dmax,
               N=counts[,.N], cbegin=cbegin, cidx=t(data.matrix(counts[,.(id1,id2)])), dist=as.array(counts[,distance]),
               eC=par$eC, log_nu=par$log_nu, log_delta=par$log_delta,
               beta_diag_centered=par$beta_diag_centered)
  capture.output(pred<-as.data.table(optimizing(stanmodels$predict_all, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)$par))
  pred=cbind(counts,pred)
  return(pred)
}

#' Fill in a large matrix with predicted values
#' @keywords internal
#' @export
#'
csnorm_predict_all_parallel = function(cs, counts, verbose=T, ncores=1) {
  nchunks=min(10*ncores,counts[,.N,by=name][,N]) #ensure at least 1 of each per parallel prediction
  counts[,chunk:=.I]
  counts[,chunk:=chunk-min(chunk),by=name]
  counts[,chunk:=as.integer(chunk/((max(chunk)+1)/nchunks)),by=name]
  registerDoParallel(cores = ncores)
  counts = foreach (i=0:(nchunks-1), .combine=rbind) %dopar%
    csnorm_predict_all(cs,counts[chunk==i],verbose=verbose)
  counts[,chunk:=NULL]
  counts
}

#' Split a chromosome into subsets for parallelization
#' @keywords internal
#'
run_split_parallel_squares = function(biases, square.size, coverage, diag.only=F) {
  minpos=biases[,min(pos)]
  maxpos=biases[,max(pos)+1]
  true.square.size=(maxpos-minpos)/round((maxpos-minpos)/square.size)
  begins=seq(minpos,maxpos-true.square.size,by=true.square.size/coverage)
  if (diag.only==T) {
    squares=SJ(begins,begins)
  } else {
    squares=CJ(begins,begins)[V2>=V1]
  }
  setnames(squares, c("begin1","begin2"))
  squares[,c("end1","end2"):=list(begin1+true.square.size,begin2+true.square.size)]
  diagsquares=squares[begin1==begin2]
  #ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2),alpha=0.1)+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=factor((begin2-begin1)/true.square.size)),alpha=0.1)+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin2-end1,xmax=end2-begin1,ymin=sqid,ymax=sqid+1),alpha=0.1)+scale_x_log10()
  return(list(squares=squares,diagsquares=diagsquares,true.square.size=true.square.size))
}

#' Initial guess for normalization
#' @keywords internal
#' @export
#' 
run_split_parallel_initial_guess = function(counts, biases, design, bf_per_kb, dmin, dmax, bf_per_decade, lambda, verbose, iter) {
  #compute column sums
  cs1=counts[,.(R=sum(contact.close+contact.down),L=sum(contact.far+contact.up)),by=c("name","id1")][,.(name,id=id1,L,R)]
  cs2=counts[,.(R=sum(contact.far+contact.down),L=sum(contact.close+contact.up)),by=c("name","id2")][,.(name,id=id2,L,R)]
  pos=biases[,.(name,id)]
  setkey(pos, name, id)
  sums=rbind(cs1,cs2)[,.(L=sum(L),R=sum(R)),keyby=c("name","id")][pos]
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)],
            XB=as.array(design[,genomic]), Krow=Krow, SD=biases[,.N], bbegin=bbegin,
            cutsitesD=biases[,pos], rejoined=biases[,rejoined],
            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
            counts_sum_left=sums[,L], counts_sum_right=sums[,R],
            lambda_nu=lambda, lambda_delta=lambda)
  op=optimizing(stanmodels$guess, data=data,
                as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
  #return initial guesses
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  return(list(eRJ=op$par$eRJ, eDE=op$par$eDE, eC=op$par$eC,
              alpha=op$par$alpha, beta_nu=op$par$beta_nu, beta_delta=op$par$beta_delta,
              log_nu=op$par$log_nu, log_delta=op$par$log_delta,
              lambda_nu=array(lambda,dim=data$Biases),
              lambda_delta=array(lambda,dim=data$Biases),
              beta_diag=beta_diag, lambda_diag=array(1,dim=Decays),
              log_decay=rep(0,counts[,.N])))
}

#' Genomic bias estimation part of parallel run
#' @keywords internal
#' 
run_split_parallel_biases = function(counts, biases, begin, end, dmin, dmax, bf_per_kb, bf_per_decade, lambda,
                                     iter, outprefix=NULL, circularize=-1) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases, begin1=begin, end1=end, fill.zeros=T, circularize=circularize)
  #initial values
  init.a=system.time(init.output <- capture.output(init.op <- run_split_parallel_initial_guess(
    counts=extracted$counts, biases=extracted$biases1, design=cs@design, bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax,
    bf_per_decade=bf_per_decade, lambda=lambda, verbose=T, iter=iter)))
  #run fit
  a=system.time(output <- capture.output(op <- csnorm_fit(extracted$biases1, extracted$counts, dmin, dmax,
                                                          bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade,
                                                          verbose = T, iter = iter, init=init.op)))
  op$runtime=a[1]+a[4]+init.a[1]+init.a[4]
  op$output=output
  op$mem=as.integer(object.size(extracted))
  if (!is.null(outprefix)) {
    save(op, file=paste0(outprefix, "_biases_op_",begin,"_",lambda,".RData"))
  }
  #report data
  center=(end-begin)/2
  ret=extracted$biases1[,.(id=id+extracted$beginrange1-1,pos,weight=(center-abs(pos-center))**2)]
  ret[,c("log_mean_DL", "log_mean_DR", "log_mean_RJ"):=list(op$par$log_mean_DL, op$par$log_mean_DR, op$par$log_mean_RJ)]
  ret[,square.begin:=begin]
  ret[,c("lambda_delta","lambda_nu","nsites","alpha"):=list(op$par$lambda_delta,op$par$lambda_nu,extracted$biases1[,.N],op$par$alpha)]
  if (!is.null(outprefix)) {
    save(ret, file=paste0(outprefix, "_biases_ret_",begin,"_",lambda,".RData"))
  }
  return(list(ret=ret, out=tail(op$output,1), runtime=op$runtime, mem=op$mem, value=op$value))
}

#' Homogenize biases estimated on multiple fits
#' @keywords internal
#' 
run_split_parallel_biases_homogenize = function(dt, lambda_nu, lambda_delta, bf_per_kb=1, iter=10000, init=0, ...) {
  Krow=round(dt[,(max(pos)-min(pos))/1000*bf_per_kb])
  data = list( Krow=Krow, S=dt[,.N], cutsites=dt[,pos],
               log_mean_RJ=dt[,log_mean_RJ], log_mean_DL=dt[,log_mean_DL], log_mean_DR=dt[,log_mean_DR],
               lambda_nu=lambda_nu, lambda_delta=lambda_delta)
  output <- capture.output(op <- optimizing(stanmodels$homogenize, data=data, as_vector=F, hessian=F,
                                            iter=iter, verbose=T, init=init, ...))
  return(op)
}

#' Decay bias estimation part of parallel run
#' @keywords internal
#' 
run_split_parallel_counts = function(counts, biases.aug, begin1, end1, begin2, end2, dmin, dmax,
                                     bf_per_decade=5, iter=100000, outprefix=NULL, circularize=-1L) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases.aug, begin1=begin1, end1=end1,
                            begin2=begin2, end2=end2, fill.zeros=T, circularize=circularize)
  #run fit
  a=system.time(output <- capture.output(op <- csnorm_fit_extradiag(
    extracted$biases1,extracted$biases2, extracted$counts,
    dmin, dmax, bf_per_decade = bf_per_decade, verbose = T, iter = iter)))
  op$runtime=a[1]+a[4]
  op$output=output
  op$mem=as.integer(object.size(extracted))
  if (!is.null(outprefix)) {
    save(op, file=paste0(outprefix, "_counts_op_",begin1,"_",begin2,".RData"))
  }
  #report data
  setkey(extracted$biases1, id, pos)
  setkey(extracted$counts, id1, pos1)
  extracted$counts[,c("log_nu1","log_delta1"):=extracted$biases1[extracted$counts,.(log_nu,log_delta)]]
  nu1off=mean(extracted$counts[,log_nu1])
  delta1off=mean(extracted$counts[,log_delta1])
  setkey(extracted$biases2, id, pos)
  setkey(extracted$counts, id2, pos2)
  extracted$counts[,c("log_nu2","log_delta2"):=extracted$biases2[extracted$counts,.(log_nu,log_delta)]]
  nu2off=mean(extracted$counts[,log_nu2])
  delta2off=mean(extracted$counts[,log_delta2])
  center=extracted$counts[,min(distance)+max(distance)]/2
  border=extracted$counts[,max(distance)-min(distance)]/2+1
  #
  setkey(extracted$counts, id1,id2,pos1,pos2,distance)
  ret=extracted$counts[,.(id1=id1+extracted$beginrange1-1,id2=id2+extracted$beginrange2-1,
                          pos1,pos2,distance,weight=(border-abs(distance-center))**2)]
  ret[,c("square.begin1","square.begin2","log_decay_close","log_decay_far","log_decay_up","log_decay_down"):=list(
    begin1, begin2,
    op$par$eC+op$par$log_decay+nu1off-delta1off+nu2off+delta2off,
    op$par$eC+op$par$log_decay+nu1off+delta1off+nu2off-delta2off,
    op$par$eC+op$par$log_decay+nu1off+delta1off+nu2off+delta2off,
    op$par$eC+op$par$log_decay+nu1off-delta1off+nu2off-delta2off)]
  if (!is.null(outprefix)) {
    save(ret, file=paste0(outprefix, "_counts_ret_",begin1,"_",begin2,".RData"))
  }
  return(list(ret=ret, out=tail(op$output,1), runtime=op$runtime, mem=op$mem))
}

#' Merge together parallel runs and estimate a single diagonal decay
#' @keywords internal
#' 
run_split_parallel_counts_decay = function(dt, counts, dmin, dmax, bf_per_decade=5, iter=10000, init=0, ...) {
  dist=dt[,(begin+end)/2]
  fij=dt[,(log_decay_close+log_decay_far+log_decay_up+log_decay_down)/4]
  ncounts=dt[,ncounts]
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=length(dist), dist=dist, fij=fij, ncounts=ncounts)
  output <- capture.output(op <- optimizing(stanmodels$fit_decay, data=data, as_vector=F,
                                            hessian=F, iter=iter, verbose=T, init=init, ...))
  return(op)
}

#' Get a proper count exposure for all estimated biases
#' @keywords internal
#' 
run_split_parallel_counts_eC = function(biases, counts, retlist, dmin, dmax, bf_per_decade=5, iter=1000) {
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               log_nu=retlist$log_nu, log_delta=retlist$log_delta, beta_diag_centered=retlist$beta_diag_centered)
  output <- capture.output(op <- optimizing(stanmodels$fit_eC, data=data, as_vector=F, hessian=F, iter=iter,
                                            verbose=T, init=0))
  return(op)
}


#' Parse intermediate outputs of killed run and try to compute the estimates nevertheless
#' @keywords internal
#' @export
#' 
run_split_parallel_recovery = function(cs, outprefix, square.size=100000, coverage=4, coverage.extradiag=1, bf_per_kb=1, bf_per_decade=5,
                                       distance_bins_per_decade=100, lambdas=c(0.01,1,100), verbose = F, iter=100000, ncores=30,
                                       homogenize=F) {
  output.binder = function(x) {
    a=sapply(names(x[[1]]), function(i) sapply(x, "[[", i,simplify=F))
    list(par=rbindlist(a[,1]), out=as.character(a[,2]), runtime=as.numeric(a[,3]))
  }
  if (file.exists(paste0(outprefix,"_ops_bias.RData"))) {
    load(paste0(outprefix,"_ops_bias.RData"))
  } else {
    ops.bias=sapply(X=Sys.glob(paste0(outprefix,"_biases_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
    ops.bias = output.binder(ops.bias)
  }
  if (file.exists(paste0(outprefix,"_ops_count.RData"))) {
    load(paste0(outprefix,"_ops_count.RData"))
  } else {
    ops.count=sapply(X=Sys.glob(paste0(outprefix,"_counts_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
    ops.count = output.binder(ops.count)
  }
  run_split_parallel(cs, design=NULL, square.size=square.size, coverage=coverage, coverage.extradiag=coverage.extradiag,
                     bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, distance_bins_per_decade=distance_bins_per_decade,
                     lambdas=lambdas, verbose=verbose, iter=iter, ncores=ncores, homogenize=homogenize, outprefix=outprefix,
                     ops.bias=ops.bias, ops.count=ops.count)
}

#' Plot bias run statistics using intermediate files
#' @keywords internal
#' @export
#' 
diagnose_biases = function(cs, outprefix, coverage=4, square.size=150000) {
  diagsquares = run_split_parallel_squares(cs@biases, square.size, coverage, diag.only=T)$diagsquares
  p0=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2),alpha=0.1)+stat_function(fun=identity)
  diagsquares[,fname:=paste0(outprefix,"_biases_op_",begin1,".RData")]
  ops = sapply(X=diagsquares[,fname], FUN=function(x){ if (file.exists(x)) {a=load(x); return(get(a[1]));} else {return(NA)}}, USE.NAMES=T, simplify=F)
  diagsquares[,runtime:=sapply(X=ops,FUN=function(x){if (class(x)=="list") {return(x$runtime)} else {return(NA)}})]
  diagsquares[,result:=as.factor(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=1) {return(tail(x$output,1)[[1]])} else {return(NA)}} else {return(NA)}}))]
  diagsquares[,nsteps:=as.integer(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=3) {return(strsplit(tail(x$output,3)[[1]],'\\s+')[[1]][2])} else {return(NA)}} else {return(NA)}}))]
  p1=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=nsteps))+stat_function(fun=identity)
  p2=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=result))+stat_function(fun=identity)
  p3=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=runtime))+stat_function(fun=identity)
  return(list(p0,p1,p2,p3,diagsquares))
}


#' Plot bias run statistics using intermediate files
#' @keywords internal
#' @export
#' 
diagnose_counts = function(cs, outprefix, coverage.extradiag=1, square.size=150000) {
  squares = run_split_parallel_squares(cs@biases, square.size, coverage.extradiag, diag.only=F)$squares
  squares[,c("id1","id2"):=list(as.integer(factor(begin1)),as.integer(factor(begin2)))]
  p0=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2))+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=is.na(nsteps)&begin2-begin1<1500000))+stat_function(fun=identity)
  squares[,fname:=paste0(outprefix,"_counts_op_",begin1,"_",begin2,".RData")]
  ops = sapply(X=squares[,fname], FUN=function(x){ if (file.exists(x)) {a=load(x); return(get(a[1]));} else {return(NA)}}, USE.NAMES=T, simplify=F)
  squares[,runtime:=sapply(X=ops,FUN=function(x){if (class(x)=="list") {return(x$runtime)} else {return(NA)}})]
  squares[,result:=as.factor(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=1) {return(tail(x$output,1)[[1]])} else {return(NA)}} else {return(NA)}}))]
  squares[,nsteps:=as.integer(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=3) {return(strsplit(tail(x$output,3)[[1]],'\\s+')[[1]][2])} else {return(NA)}} else {return(NA)}}))]
  p1=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=nsteps))+stat_function(fun=identity)
  p2=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=result))+stat_function(fun=identity)
  p3=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=runtime))+stat_function(fun=identity)
  return(list(p0,p1,p2,p3,squares))
}

#' Cut-site normalization (parallelized)
#' 
#' The dataset is splitted into subsets which are run in parallel and they are 
#' then stitched together. First, genomic biases nu and delta are estimated by 
#' cutting the genome in overlapping segments. Mean estimates are then computed 
#' for nu and delta. These estimates can then be homogenized to have a constant 
#' stiffness. Second, squares covering all the dataset are extracted in order to
#' compute the diagonal decay given nu and delta. Again, estimates are computed 
#' for the diagonal decay. Finally, the count exposure is recomputed for the 
#' final estimates of nu, delta and decay.
#' 
#' @section Note: This call will spawn parallel processes using foreach and
#'   doParallel. The number of parallel processes is set by ncores.
#'   
#' @param cs CSnorm object as returned by \code{\link{merge_cs_norm_datasets}}
#' @param square.size positive integer. Size of the subset of data (in base 
#'   pairs) to normalize independently. If too large, optimization fails to 
#'   converge and takes too long. If too small, estimates are heavily biased.
#' @param coverage positive integer. To reduce border effects, squares overlap 
#'   and estimates are computed using a weighted average that favors the center 
#'   of the square. Coverage indicates how many squares cover any portion of the
#'   genome. This parameter is for the estimation of the genomic biases.
#' @param coverage.extradiag positive integer. Same as previous, but in the 
#'   decay step.
#' @param bf_per_kb positive numeric. Number of cubic spline basis functions per
#'   kilobase, for genomic bias estimates. Small values make the optimization 
#'   easy, but makes the genomic biases stiffer.
#' @param bf_per_decade positive numeric. Number of cubic spline basis functions
#'   per distance decade (in bases), for diagonal decay. Default parameter 
#'   should suffice.
#' @param distance_bins_per_decade positive integer. How many bins per decade 
#'   should be used to discretize and re-fit the diagonal decay. Default 
#'   parameter should suffice. Must be much larger than bf_per_decade: ideally 
#'   in that bin any basis function should be approximately constant.
#' @param lambdas vector of positive values. Length scales to try out as initial
#'   conditions.
#' @param verbose boolean. Show output of different steps.
#' @param iter positive integer. Number of optimization steps for each stan 
#'   model call.
#' @param ncores positive integer. Number of cores to parallelize on.
#' @param homogenize boolean. Should the biases be homogenized for stiffness? 
#'   Default is FALSE
#' @param outprefix character. If not NULL, prefix used to write intermediate 
#'   output files. Diagnostics only.
#' @param circularize integer. Set this to the size of the chromosome if it is 
#'   circular, otherwise leave as-is (default is -1)
#' @param ops.bias,ops.count if not NULL, skip the genomic (resp. decay) bias 
#'   estimation step, and use these intermediate data tables instead.
#'   
#' @return A list containing: \enumerate{ \item par: the optimized parameters 
#'   \item out.bias: the last line of the stan output for the bias estimation 
#'   \item runtime.bias: the runtimes of these optimizations \item out.count:
#'   the last line of the stan output for the decay estimation \item
#'   runtime.count: the runtimes of these optimizations }
#' @export
#' 
#' @examples
run_split_parallel = function(cs, design=NULL, square.size=100000, coverage=4, coverage.extradiag=1, bf_per_kb=1, bf_per_decade=5,
                              distance_bins_per_decade=100, lambdas=c(0.01,1,100), verbose = F, iter=100000, ncores=30, homogenize=F,
                              outprefix=NULL, ops.bias=NULL, ops.count=NULL) {
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]/2))
  cs@settings = c(cs@settings, list(square.size=square.size, coverage=coverage, coverage.extradiag=coverage.extradiag,
                                    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, distance_bins_per_decade=distance_bins_per_decade,
                                    iter=iter, homogenize=homogenize, outprefix=outprefix))
  stopifnot(is.null(design))
  ### build squares
  if (verbose==T) cat("*** build squares\n")
  retval.squares = run_split_parallel_squares(cs@biases, square.size, coverage, diag.only=T)
  diagsquares=retval.squares$diagsquares
  true.square.size=retval.squares$true.square.size
  retval.squares = run_split_parallel_squares(cs@biases, square.size, coverage.extradiag, diag.only=F)
  stopifnot(true.square.size==retval.squares$true.square.size)
  squares=retval.squares$squares
  if (verbose==T) {
    cat("Run in parallel\n")
    cat("   square size: ",true.square.size,"\n")
    cat("   diagonal coverage: ",coverage,"\n")
    cat("   extradiagonal coverage: ",coverage.extradiag,"\n")
    cat("   ",squares[,.N], " squares\n")
    cat("   ",diagsquares[,.N], " on the diagonal\n")
  }
  ### fit genomic biases using squares close to diagonal
  if (verbose==T) cat("*** fit genomic biases\n")
  dmin=0.99
  if (cs@settings$circularize>0) {
    dmax=cs@settings$circularize/2+0.01
  } else {
    dmax=cs@biases[,max(pos)-min(pos)]+0.01
  }
  cs@settings$dmin=dmin
  cs@settings$dmax=dmax
  registerDoParallel(cores=ncores)
  output.binder = function(x) {
    a=sapply(names(x[[1]]), function(i) sapply(x, "[[", i,simplify=F))
    list(par=rbindlist(a[,1]), out=as.character(a[,2]), runtime=as.numeric(a[,3]), mem=as.integer(a[,4]))
  }
  if (is.null(ops.bias)) {
    ops.bias = foreach (begin=diagsquares[,begin1], end=diagsquares[,end1], .packages=c("data.table","rstan","csnorm")) %:%
      foreach (lambda=lambdas, .combine=function(x,y){if (x$value<y$value) {return(y)} else {return(x)}}) %dopar% 
      run_split_parallel_biases(cs@counts, cs@biases, begin, end, dmin, dmax, bf_per_kb = bf_per_kb,
                                bf_per_decade = bf_per_decade, lambda=lambda, iter = iter, outprefix=outprefix,
                                circularize=cs@settings$circularize)
    ops.bias = output.binder(ops.bias)
    if (!is.null(outprefix)) {save(ops.bias, file=paste0(outprefix, "_ops_bias.RData"))}
  }
  #ggplot(ops.bias$par)+geom_line(aes(pos,log_mean_RJ,colour=factor(square.begin)))+guides(colour=F)
  ### reconstruct bias estimates: eRJ eDE log_nu and log_delta
  if (verbose==T) cat("*** reconstruct genomic biases\n")
  info=ops.bias$par[,.SD[1],by=square.begin,.SDcols=c("lambda_nu","lambda_delta","nsites","alpha")]
  #ggplot(info)+geom_line(aes(square.begin,lambda_nu))+scale_y_log10()+geom_hline(yintercept=info[,exp(median(log(lambda_nu)))])
  means=ops.bias$par[,.(log_mean_RJ=weighted.mean(log_mean_RJ, weight),
                        log_mean_DL=weighted.mean(log_mean_DL, weight),
                        log_mean_DR=weighted.mean(log_mean_DR, weight),
                        ncounts=.N), keyby=c("id","pos")]
  if (homogenize==T) {
    if (verbose==T) cat("*** homogenize genomic biases\n")
    op=run_split_parallel_biases_homogenize(means, lambda_nu=info[,exp(mean(log(lambda_nu)))], lambda_delta=info[,exp(mean(log(lambda_delta)))],
                                            bf_per_kb=bf_per_kb, iter=iter)
    #ggplot(means)+geom_line(aes(pos,log_mean_RJ,colour="ori"))+
    #  geom_line(data=data.table(pos=means[,pos],rj=op$par$log_nu+op$par$eRJ),aes(pos,rj,colour="avg"))
    #ggplot(means)+geom_line(aes(pos,log_mean_DL,colour="ori"))+
    #  geom_line(data=data.table(pos=means[,pos],rj=op$par$log_nu+op$par$eDE+op$par$log_delta),aes(pos,rj,colour="avg"))
    retlist=op$par[c("eRJ","eDE","log_nu","log_delta")]
    setkey(cs@biases, id, pos)
    biases.aug=copy(cs@biases)
    biases.aug[,log_nu:=retlist$log_nu]
    biases.aug[,log_delta:=retlist$log_delta]
  } else {
    setkey(cs@biases, id, pos)
    retlist=list(eRJ=ops.bias$par[,mean(log_mean_RJ)], eDE=ops.bias$par[,mean(log_mean_DL+log_mean_DR)/2])
    biases.aug=means[cs@biases]
    stopifnot(!any(is.na(biases.aug)))
    biases.aug[,log_nu:=log_mean_RJ-retlist$eRJ]
    biases.aug[,log_delta:=(log_mean_DL-log_mean_DR)/2]
    retlist$log_nu=biases.aug[,log_nu]
    retlist$log_delta=biases.aug[,log_delta]
  }
  #ggplot(melt(biases.aug[,.(pos,log_nu,log_delta)], id.vars="pos"))+geom_line(aes(pos,value,colour=variable))
  #ggplot(melt(biases.aug[,.(pos,nu=exp(log_nu),delta=exp(log_delta))], id.vars="pos"))+geom_line(aes(pos,value,colour=variable))
  #
  ### fit remaining data
  if (verbose==T) cat("*** fit diagonal decay\n")
  registerDoParallel(cores=ncores)
  if (is.null(ops.count)) {
    ops.count = foreach (begin1=squares[,begin1], end1=squares[,end1], begin2=squares[,begin2], end2=squares[,end2],
                         .packages=c("data.table","rstan")) %dopar% 
      run_split_parallel_counts(cs@counts, biases.aug, begin1, end1, begin2, end2, dmin, dmax,
                                bf_per_decade = bf_per_decade, iter = iter, outprefix=outprefix,
                                circularize=cs@settings$circularize)
    ops.count=output.binder(ops.count)
    if (!is.null(outprefix)) {save(ops.count, file=paste0(outprefix, "_ops_count.RData"))}
  }
  ### reconstruct count estimates: log_decay and beta_diag_centered
  if (verbose==T) cat("*** reconstruct diagonal decay\n")
  #qplot(ops.count[,id],binwidth=1)
  #ops.count[,cat:=paste(square.begin1,square.begin2)]
  #ops.count[,dbin:=factor(square.begin2-square.begin1)]
  #ggplot(ops.count)+geom_line(aes(distance,weight,colour=cat))+guides(colour=F)+facet_grid(dbin~.)
  #ggplot(ops.count[weight>1e9])+geom_line(aes(distance,log_decay_far,colour=cat))+guides(colour=F)+facet_grid(.~dbin)
  #ggplot(ops.count[square.begin1==square.begin2&square.begin1==3189.0])+geom_line(aes(distance,log_decay_up))
  #
  stepsz=1/distance_bins_per_decade
  #if (circularize>0) {
  #  dmax=biases.aug[,max(pmin(pos,circularize-pos+1))]
  #  dmin=biases.aug[,min(pmin(pos,circularize-pos+1))]
  #  dbins=10**seq(0,log10(dmax-dmin)+stepsz,stepsz)
  #} else {
  dbins=10**seq(log10(dmin),log10(dmax)+stepsz,stepsz)
  #}
  ops.count$par[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #ops.count$par[,grid:=factor(square.begin2-square.begin1)]
  #ggplot(ops.count$par[sample(.N,100000)])+geom_line(aes(distance,log_decay_far,group=square.begin1-square.begin2,colour=bdist))+guides(colour=F)+facet_grid(.~grid)
  test=ops.count$par[,.(log_decay_close=weighted.mean(log_decay_close, weight),
                        log_decay_far=weighted.mean(log_decay_far, weight),
                        log_decay_up=weighted.mean(log_decay_up, weight),
                        log_decay_down=weighted.mean(log_decay_down, weight),
                        ncounts=.N), keyby=bdist] #could only take one decay, and even weights since eC is estimated later
  bin.begin=test[,bdist]
  bin.end=test[,bdist]
  levels(bin.begin) <- tstrsplit(as.character(levels(bin.begin)), "[[,]")[2][[1]]
  levels(bin.end) <- tstrsplit(as.character(levels(bin.end)), "[],)]")[2][[1]]
  test[,begin:=as.integer(as.character(bin.begin))]
  test[,end:=as.integer(as.character(bin.end))]
  op=run_split_parallel_counts_decay(test, counts, dmin, dmax, bf_per_decade = bf_per_decade, iter = iter)
  test[,log_decay:=op$par$log_decay+op$par$eC]
  retlist$decay=test[,.(dist=(begin+end)/2,decay=exp(log_decay))]
  retlist$beta_diag_centered=op$par$beta_diag_centered
  ### reconstruct count estimates: eC
  if (verbose==T) cat("*** reconstruct count exposure\n")
  op=run_split_parallel_counts_eC(cs@biases, cs@counts[sample(.N,min(.N,100000))], retlist, dmin, dmax,
                                  bf_per_decade=bf_per_decade)
  retlist$eC=op$par$eC
  retlist$alpha=op$par$alpha
  #ggplot(test)+geom_point(aes(bdist, log_decay))
  #ggplot(melt(test, id.vars=c("bdist","begin","end")))+geom_point(aes(bdist,value,colour=variable))
  if (verbose==T) cat("*** done\n")
  cs@par=retlist
  cs@diagnostics=list(out.bias=ops.bias$out, out.count=ops.count$out, runtime.count=ops.count$runtime, runtime.bias=ops.bias$runtime)
  return(cs)
}

#' Run approximate gibbs sampler on with a single starting condition
#' @inheritParams run_simplified
#' @param fit.decay,fit.genomic boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @keywords internal
#' @export
#' 
run_simplified_gibbs = function(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, lambda=1,
                                ngibbs = 3, iter=100000, fit.decay=T, fit.genomic=T, verbose=T, ncounts=100000) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf, groups=groups,
                                    lambda=lambda, iter=iter, ngibbs=ngibbs))
  #fill counts matrix and take subset
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize)
  if (cs@counts[,.N] > ncounts*2+1) {
    subcounts = cs@counts[,.SD[(1:ncounts+as.integer(.N/2))],by=name]
  } else {
    subcounts=cs@counts
  }
  if (subcounts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2) stop("dataset too sparse, please increase ncounts")
  #report min/max distance
  dmin=0.99
  if (cs@settings$circularize>0) {
    dmax=cs@settings$circularize/2+0.01
  } else {
    dmax=cs@biases[,max(pos)-min(pos)]+0.01
  }
  cs@settings$dmin=dmin
  cs@settings$dmax=dmax
  #initial guess
  if (verbose==T) cat("Initial guess\n")
  init.a = system.time(init.output <- capture.output(init.op <- csnorm:::csnorm_simplified_guess(
    biases = cs@biases, counts = cs@counts, design = cs@design, lambda=lambda, dmin=dmin, dmax=dmax,
    groups = groups, bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, iter = iter)))
  cs@diagnostics=list(out.init=init.output, runtime.init=init.a[1]+init.a[4])
  #abort silently if initial guess went wrong
  if (length(grep("Line search failed",tail(init.output,1)))>0) {
    init.op$par$value=.Machine$double.xmax
    cs@par=init.op$par
    return(cs)
  }
  op=init.op
  #make sure beta_diag is strictly increasing
  for (d in 2:length(op$par$beta_diag)) {
    if (abs(op$par$beta_diag[d]-op$par$beta_diag[d-1])<10*.Machine$double.eps) {
      op$par$beta_diag[d:length(op$par$beta_diag)]=op$par$beta_diag[d:length(op$par$beta_diag)]+10*.Machine$double.eps
    }
  }
  #gibbs sampling
  for (i in 1:ngibbs) {
    #fit diagonal decay given nu and delta
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay\n")
      a=system.time(output <- capture.output(op.diag <- csnorm:::csnorm_simplified_decay(
        biases = cs@biases, counts = cs@counts, design=cs@design,
        log_nu = op$par$log_nu, log_delta = op$par$log_delta,
        dmin = dmin, dmax = dmax, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, groups = groups,
        iter=iter, init=op$par, dispersion=op$par$alpha)))
      op=list(value=op.diag$value, par=c(op.diag$par[c("eC","beta_diag","beta_diag_centered",
                                                       "lambda_diag","log_decay","decay")],
                                         op$par[c("eRJ","eDE","beta_nu","beta_delta", "alpha",
                                                  "lambda_nu","lambda_delta","log_nu","log_delta")]))
      cs@diagnostics[[paste0("out.decay",i)]]=output
      cs@diagnostics[[paste0("runtime.decay",i)]]=a[1]+a[4]
    }
    #fit nu and delta given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic\n")
      a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_simplified_genomic(
        biases = cs@biases, counts = cs@counts, design = cs@design,
        log_decay = op$par$log_decay, log_nu = op$par$log_nu, log_delta = op$par$log_delta,
        groups = groups, bf_per_kb = bf_per_kb, iter = iter, init=op$par, dispersion=op$par$alpha)))
      op=list(value=op.gen$value, par=c(op$par[c("beta_diag","beta_diag_centered","lambda_diag","log_decay","decay", "alpha")],
                                        op.gen$par[c("eC","eRJ","eDE","beta_nu","beta_delta",
                                                     "lambda_nu","lambda_delta","log_nu","log_delta")]))
      cs@diagnostics[[paste0("out.bias",i)]]=output
      cs@diagnostics[[paste0("runtime.bias",i)]]=a[1]+a[4]
    }
    #make sure beta_diag is strictly increasing
    for (d in 2:length(op$par$beta_diag)) {
      if (abs(op$par$beta_diag[d]-op$par$beta_diag[d-1])<10*.Machine$double.eps) {
        op$par$beta_diag[d:length(op$par$beta_diag)]=op$par$beta_diag[d:length(op$par$beta_diag)]+10*.Machine$double.eps
      }
    }
    #fit exposures and dispersion
    a=system.time(output <- capture.output(op.disp <- csnorm:::csnorm_simplified_dispersion(
      biases = cs@biases, counts = subcounts, design = cs@design,
      dmin=dmin, dmax=dmax, beta_nu=op$par$beta_nu, beta_delta=op$par$beta_delta, beta_diag=op$par$beta_diag,
      lambda_nu=op$par$lambda_nu, lambda_delta=op$par$lambda_delta, lambda_diag=op$par$lambda_diag,
      bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=iter, init=op$par)))
    op=list(value=op.disp$value, par=c(op$par[c("beta_diag","beta_diag_centered","lambda_diag","log_decay","decay",
                                                "beta_nu","beta_delta", "lambda_nu","lambda_delta","log_nu","log_delta")],
                                      op.disp$par[c("eC","eRJ","eDE","alpha")]))
    cs@diagnostics[[paste0("out.disp",i)]]=output
    cs@diagnostics[[paste0("runtime.disp",i)]]=a[1]+a[4]
  }
  if (verbose==T) cat("Done\n")
  op$par$runtime=sum(as.numeric(cs@diagnostics[grep("runtime",names(cs@diagnostics))]))
  op$par$output=output
  init.op$par$runtime=init.a[1]+init.a[4]
  init.op$par$output=init.output
  op$par$init=init.op$par
  op$par$value=op$value
  cs@par=op$par
  return(cs)
}

#' Cut-site normalization (simplified gibbs sampling)
#' 
#' Alternates two approximations to the exact model, fitting the diagonal decay and nu/delta.
#' 
#' @param cs CSnorm object as returned by \code{\link{merge_cs_norm_datasets}}
#' @param bf_per_kb positive numeric. Number of cubic spline basis functions per
#'   kilobase, for genomic bias estimates. Small values make the optimization 
#'   easy, but makes the genomic biases stiffer.
#' @param bf_per_decade positive numeric. Number of cubic spline basis functions
#'   per distance decade (in bases), for diagonal decay. Default parameter 
#'   should suffice.
#' @param bins_per_bf positive integer. Number of distance bins to split basis
#'   functions into. Must be sufficiently small so that the diagonal decay is
#'   approximately constant in that bin.
#' @param groups positive integer. Number of quantile groups to keep in each
#'   distance bin.
#' @param lambdas positive numeric. Length scales to try out as initial condition.
#' @param ngibbs positive integer. Number of gibbs sampling iterations.
#' @param iter positive integer. Number of optimization steps for each stan
#'   optimization call.
#' @param ncores positive integer. Number of cores to parallelize on.
#' @param verbose Display progress if TRUE
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
run_simplified = function(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=c(0.1,1,10),
                          ngibbs = 3, iter=100000, ncores=1, verbose=T) {
  cs@binned=list() #erase old binned datasets if available
  registerDoParallel(cores=ncores)
  cs = foreach (lambda=lambdas, .combine=function(x,y){if (x@par$value[1]<y@par$value[1]){return(y)}else{return(x)}}) %dopar%
    run_simplified_gibbs(cs, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade,
                         bins_per_bf=bins_per_bf, groups=groups, lambda=lambda, ngibbs = ngibbs,
                         iter=iter, fit.decay=T, fit.genomic=T, verbose=verbose)
  return(cs)
}


#' Run approximate gibbs sampler on with a single starting condition
#' @inheritParams run_simplified
#' @param fit.decay,fit.genomic boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @keywords internal
#' @export
#' 
run_gauss_gibbs = function(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, lambda=1,
                                ngibbs = 3, iter=100000, fit.decay=T, fit.genomic=T, verbose=T, ncounts=100000) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
                                    lambda=lambda, iter=iter, ngibbs=ngibbs))
  #fill counts matrix and take subset
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize)
  if (cs@counts[,.N] > ncounts*2+1) {
    subcounts = cs@counts[,.SD[(1:ncounts+as.integer(.N/2))],by=name]
  } else {
    subcounts=cs@counts
  }
  if (subcounts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2) stop("dataset too sparse, please increase ncounts")
  #report min/max distance
  dmin=0.99
  if (cs@settings$circularize>0) {
    dmax=cs@settings$circularize/2+0.01
  } else {
    dmax=cs@biases[,max(pos)-min(pos)]+0.01
  }
  cs@settings$dmin=dmin
  cs@settings$dmax=dmax
  #initial guess
  if (verbose==T) cat("Initial guess\n")
  init.a = system.time(init.output <- capture.output(init.op <- csnorm:::csnorm_gauss_guess(
    biases = cs@biases, counts = cs@counts, design = cs@design, lambda=lambda, dmin=dmin, dmax=dmax,
    bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, iter = iter)))
  cs@diagnostics=list(out.init=init.output, runtime.init=init.a[1]+init.a[4])
  #abort silently if initial guess went wrong
  if (length(grep("Line search failed",tail(init.output,1)))>0) {
    init.op$par$value=.Machine$double.xmax
    cs@par=init.op$par
    return(cs)
  }
  op=init.op
  #make sure beta_diag is strictly increasing
  for (d in 2:length(op$par$beta_diag)) {
    if (abs(op$par$beta_diag[d]-op$par$beta_diag[d-1])<10*.Machine$double.eps) {
      op$par$beta_diag[d:length(op$par$beta_diag)]=op$par$beta_diag[d:length(op$par$beta_diag)]+10*.Machine$double.eps
    }
  }
  #gibbs sampling
  for (i in 1:ngibbs) {
    #fit diagonal decay given nu and delta
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay\n")
      a=system.time(output <- capture.output(op.diag <- csnorm:::csnorm_gauss_decay(
        biases = cs@biases, counts = cs@counts, design=cs@design,
        log_nu = op$par$log_nu, log_delta = op$par$log_delta,
        dmin = dmin, dmax = dmax, dispersion=op$par$alpha, lambda_diag=op$par$lambda_diag,
        bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf,
        iter=iter, init=op$par)))
      op=list(value=op.diag$value, par=c(op.diag$par[c("eC","beta_diag","beta_diag_centered",
                                                       "log_decay","decay")],
                                         op$par[c("eRJ","eDE","beta_nu","beta_delta", "alpha","lambda_diag",
                                                  "lambda_nu","lambda_delta","log_nu","log_delta")]))
      cs@diagnostics[[paste0("out.decay",i)]]=output
      cs@diagnostics[[paste0("runtime.decay",i)]]=a[1]+a[4]
    }
    #fit nu and delta given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic\n")
      a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_gauss_genomic(
        biases = cs@biases, counts = cs@counts, design = cs@design,
        log_decay = op$par$log_decay, log_nu = op$par$log_nu, log_delta = op$par$log_delta,
        dispersion=op$par$alpha, lambda_nu=op$par$lambda_nu, lambda_delta=op$par$lambda_delta,
        bf_per_kb = bf_per_kb, iter = iter, init=op$par)))
      op=list(value=op.gen$value, par=c(op.gen$par[c("eC","eRJ","eDE","beta_nu","beta_delta",
                                                     "log_nu","log_delta")],
                                        op$par[c("beta_diag","beta_diag_centered","lambda_diag",
                                                 "lambda_nu","lambda_delta","log_decay","decay", "alpha")]))
      cs@diagnostics[[paste0("out.bias",i)]]=output
      cs@diagnostics[[paste0("runtime.bias",i)]]=a[1]+a[4]
    }
    #make sure beta_diag is strictly increasing
    for (d in 2:length(op$par$beta_diag)) {
      if (abs(op$par$beta_diag[d]-op$par$beta_diag[d-1])<10*.Machine$double.eps) {
        op$par$beta_diag[d:length(op$par$beta_diag)]=op$par$beta_diag[d:length(op$par$beta_diag)]+10*.Machine$double.eps
      }
    }
    #fit exposures and dispersion
    a=system.time(output <- capture.output(op.disp <- csnorm:::csnorm_gauss_dispersion(
      biases = cs@biases, counts = subcounts, design = cs@design,
      dmin=dmin, dmax=dmax, beta_nu=op$par$beta_nu, beta_delta=op$par$beta_delta, beta_diag=op$par$beta_diag,
      bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=iter, init=op$par)))
    op=list(value=op.disp$value, par=c(op.disp$par[c("eC","eRJ","eDE","alpha",
                                                     "lambda_diag","lambda_nu","lambda_delta")],
                                       op$par[c("beta_diag","beta_diag_centered", "log_decay","decay",
                                                "beta_nu","beta_delta", "log_nu","log_delta")]))
    cs@diagnostics[[paste0("out.disp",i)]]=output
    cs@diagnostics[[paste0("runtime.disp",i)]]=a[1]+a[4]
  }
  if (verbose==T) cat("Done\n")
  op$par$runtime=sum(as.numeric(cs@diagnostics[grep("runtime",names(cs@diagnostics))]))
  op$par$output=output
  init.op$par$runtime=init.a[1]+init.a[4]
  init.op$par$output=init.output
  op$par$init=init.op$par
  op$par$value=op$value
  cs@par=op$par
  return(cs)
}

#' Cut-site normalization (simplified gibbs sampling)
#' 
#' Alternates two approximations to the exact model, fitting the diagonal decay and nu/delta.
#' 
#' @param cs CSnorm object as returned by \code{\link{merge_cs_norm_datasets}}
#' @param bf_per_kb positive numeric. Number of cubic spline basis functions per
#'   kilobase, for genomic bias estimates. Small values make the optimization 
#'   easy, but makes the genomic biases stiffer.
#' @param bf_per_decade positive numeric. Number of cubic spline basis functions
#'   per distance decade (in bases), for diagonal decay. Default parameter 
#'   should suffice.
#' @param bins_per_bf positive integer. Number of distance bins to split basis
#'   functions into. Must be sufficiently small so that the diagonal decay is
#'   approximately constant in that bin.
#' @param lambdas positive numeric. Length scales to try out as initial condition.
#' @param ngibbs positive integer. Number of gibbs sampling iterations.
#' @param iter positive integer. Number of optimization steps for each stan
#'   optimization call.
#' @param ncores positive integer. Number of cores to parallelize on.
#' @param verbose Display progress if TRUE
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
run_gauss = function(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, lambdas=c(0.1,1,10),
                          ngibbs = 3, iter=100000, ncores=1, verbose=T) {
  cs@binned=list() #erase old binned datasets if available
  registerDoParallel(cores=ncores)
  cs = foreach (lambda=lambdas, .combine=function(x,y){if (x@par$value[1]<y@par$value[1]){return(y)}else{return(x)}}) %dopar%
    run_gauss_gibbs(cs, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade,
                         bins_per_bf=bins_per_bf, lambda=lambda, ngibbs = ngibbs,
                         iter=iter, fit.decay=T, fit.genomic=T, verbose=verbose)
  return(cs)
}


#' Run exact model on a single cpu
#' @inheritParams run_exact
#' @keywords internal
#' @export
#' 
run_serial = function(cs, bf_per_kb=1, bf_per_decade=5, lambda=1, iter=100000, subsampling.pc=100) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, lambda=lambda, iter=iter))
  #fill counts matrix
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize)
  #report min/max distance
  dmin=0.99
  if (cs@settings$circularize>0) {
    dmax=cs@settings$circularize/2+0.01
  } else {
    dmax=cs@biases[,max(pos)-min(pos)]+0.01
  }
  cs@settings$dmin=dmin
  cs@settings$dmax=dmax
  setkey(cs@biases,name,id,pos)
  setkey(cs@counts,name,id1,pos1,id2,pos2)
  #initial guess
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases, design=cs@design,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, lambda=lambda, verbose=T, iter=iter)))
  cs@diagnostics=list(out.init=init.output, runtime.init=init.a[1]+init.a[4])
  #main optimization, subsampled
  counts.sub=cs@counts[sample(.N,round(subsampling.pc/100*.N))]
  setkeyv(counts.sub,key(cs@counts))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = counts.sub, design=cs@design, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=iter, verbose = T,
    init=init.op, weight=counts.sub[,.N,by=name]$N/cs@counts[,.N,by=name]$N)))
  cs@diagnostics=list(out=output, runtime=a[1]+a[4])
  #report statistics
  op$par$init=init.op
  op$par$value=op$value
  if (subsampling.pc<100) op$par$counts.sub=counts.sub
  cs@par=op$par
  cs
}



#' Cut-site normalization (exact model)
#' 
#' Will run the exact model of normalization (on one cpu for each lambda 
#' provided) and returns the most likely model and predicted quantities. Useful
#' for comparison purposes. If you don't know what to use, try 
#' \code{\link{run_simplified}}.
#' 
#' @inheritParams run_simplified
#' @param subsampling.pc numeric. Percentage of the data used to do the calculations (default 100).
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
run_exact = function(cs, bf_per_kb=1, bf_per_decade=5, lambdas=c(0.1,1,10), ncores=1, iter=100000, subsampling.pc=100) {
  cs@binned=list() #erase old binned datasets if available
  registerDoParallel(cores=ncores)
  cs = foreach (lambda=lambdas, .combine=function(x,y){if (x@par$value<y@par$value){return(y)}else{return(x)}}) %dopar%
    run_serial(cs, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, lambda=lambda, iter=iter, subsampling.pc=subsampling.pc)
  return(cs)
}

#' Verify model fit by computing posterior predictive quantities
#'
#' @param cs CSnorm object, normalized.
#' @param genomic.groups How many groups for nu and delta
#' @param decay.groups How many groups for diagonal decay
#' @param npoints Number of points per group to take (not guaranteed)
#'
#' @return
#' @export
#'
#' @examples
check_fit = function(cs, genomic.groups=5, decay.groups=5, npoints=10) {
  if (length(cs@par)==0) stop("Must normalize the datasets first")
  #build bins
  dbins=c(0,10**seq(3,log10(cs@settings$dmax),length.out=decay.groups))
  gbins=cut2(c(cs@par$log_nu,cs@par$log_delta),g=genomic.groups,onlycuts=T)
  #build counts matrix 
  biases=copy(cs@biases[,.(name,id,pos)])
  biases[,c("log_nu","log_delta"):=list(cs@par$log_nu,cs@par$log_delta)]
  biases[,nubin:=cut(log_nu, gbins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  biases[,deltabin:=cut(log_delta, gbins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  biases[,c("log_nu","log_delta"):=list(NULL,NULL)]
  biases=biases[,.SD[sample(.N,min(.N,npoints))],by=c("name","nubin","deltabin")]
  counts=cs@counts[id1%in%biases[,id]&id2%in%biases[,id]]
  counts=fill_zeros(counts,biases,circularize=cs@settings$circularize)
  #filter by distance
  counts[,dbin:=cut(distance, dbins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  counts=counts[,.SD[sample(.N,min(.N,npoints))],by=c("name","dbin")]
  setkey(counts, name, id1, id2, pos1, pos2)
  #predict values
  counts=csnorm_predict_all(cs, counts, verbose=F)
  counts=rbind(counts[,.(name,id1,id2,dbin,count=contact.close,mean=exp(log_mean_cclose))],
               counts[,.(name,id1,id2,dbin,count=contact.far,mean=exp(log_mean_cfar))],
               counts[,.(name,id1,id2,dbin,count=contact.up,mean=exp(log_mean_cup))],
               counts[,.(name,id1,id2,dbin,count=contact.down,mean=exp(log_mean_cdown))])
  biases[,pos:=NULL]
  counts=merge(counts,biases,by.x=c("name","id1"),by.y=c("name","id"))
  counts=merge(counts,biases,by.x=c("name","id2"),by.y=c("name","id"),suffixes=c("1","2"))
  #compute p-values
  counts[count>=mean,pval:=pnbinom(count,size=cs@par$alpha,mu=mean,lower.tail=F)]
  counts[count<mean,pval:=pnbinom(count,size=cs@par$alpha,mu=mean,lower.tail=T)]
  counts[,sd:=sqrt(mean+mean**2/cs@par$alpha)]
  #graph p-values
  p.all=ggplot(counts)+geom_histogram(aes(pval))+facet_wrap(~name)+xlab("model p-value")+ylab("frequency")
  p.decay=ggplot(counts)+geom_jitter(aes(1,pval))+facet_grid(name~dbin)+ylab("model p-value")+xlab("diagonal decay bin")
  p.nu=ggplot(counts)+geom_jitter(aes(1,pval))+facet_grid(name~nubin1)+ylab("model p-value")+xlab("nu bin")
  p.delta=ggplot(counts)+geom_jitter(aes(1,pval))+facet_grid(name~deltabin1)+ylab("model p-value")+xlab("delta bin")
  return(list(all=p.all,decay=p.decay,nu=p.nu,delta=p.delta,counts=counts))
}



