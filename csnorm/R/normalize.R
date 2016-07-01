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
fill_zeros = function(counts,biases,biases2=NULL) {
  if (is.null(biases2)) biases2=biases
  newcounts=CJ(biases[,paste(id,pos)],biases2[,paste(id,pos)])
  newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
  newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
  newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
  newcounts=newcounts[pos1<pos2]
  setkey(newcounts, id1, id2, pos1, pos2)
  setkey(counts, id1, id2, pos1, pos2)
  newcounts=counts[newcounts]
  newcounts[is.na(contact.close),contact.close:=0]
  newcounts[is.na(contact.far),contact.far:=0]
  newcounts[is.na(contact.up),contact.up:=0]
  newcounts[is.na(contact.down),contact.down:=0]
  return(newcounts)
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
  if (fill.zeros==T) {
    counts.local=fill_zeros(counts=counts.local,biases=biases.1,biases2=biases.2)
    if (circularize>0) {
      counts.local[,distance:=pmin(abs(pos2-pos1), circularize+1-abs(pos2-pos1))]
    } else {
      counts.local[,distance:=abs(pos2-pos1)]
    }
  }
  return(list(counts=counts.local,biases1=biases.1,biases2=biases.2,
              beginrange1=beginrange1, endrange1=endrange1, beginrange2=beginrange2, endrange2=endrange2))
}

#' Single-cpu fitting
#' @keywords internal
#' 
csnorm_fit = function(model, biases, counts, dmin, dmax, bf_per_kb=1, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far], counts_up=counts[,contact.up], counts_down=counts[,contact.down])
  if (verbose==T) {
    message("CS norm: fit")
    message("Krow        : ", Krow)
    message("Kdiag       : ", Kdiag)
    message("Biases      : ", biases[,.N])
    message("Counts      : ", 4*counts[,.N])
    message("% zeros     : ", (counts[contact.close==0,.N]+counts[contact.far==0,.N]+counts[contact.up==0,.N]+counts[contact.down==0,.N])/counts[,.N]*100)
  }
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

#' Single-cpu fitting, only diagonal decay
#' @keywords internal
#' 
csnorm_fit_extradiag = function(model, biases1, biases2, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( S1=biases1[,.N], S2=biases2[,.N], Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far], counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               log_nu1=biases1[,log_nu], log_nu2=biases2[,log_nu], log_delta1=biases1[,log_delta], log_delta2=biases2[,log_delta])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

#' Predict expected values for each count given optimized model parameters
#'
#' @param cs CSnorm object
#' @param verbose
#' @param ncores
#'
#' @return a stan optimization output with the predictions.
#' @keywords internal
#' @export
#'
#' @examples
csnorm_predict_all = function(cs, verbose=T, ncores=1) {
  biases=cs@biases
  par=cs@par
  dmin=cs@settings$dmin
  dmax=cs@settings$dmax
  Kdiag=round((log10(dmax)-log10(dmin))*cs@settings$bf_per_decade)
  registerDoParallel(cores=ncores)
  ncounts=cs@counts[,.N]
  stepsize=ncores*10
  output.binder = function(x) {
    a=sapply(names(x[[1]]), function(i) sapply(x, "[[", i,simplify=F))
    list(par=rbindlist(a[,1]), out=as.character(a[,2]), runtime=as.numeric(a[,3]))
  }
  pred = foreach (i=seq(1,ncounts+stepsize,stepsize), .packages=c("data.table","rstan"), .combine="rbind") %dopar% {
    counts=cs@counts[i:min(.N,i+stepsize-1)]
    cclose=counts[,.(id1,id2,distance,count=contact.close)]
    cfar=counts[,.(id1,id2,distance,count=contact.far)]
    cup=counts[,.(id1,id2,distance,count=contact.up)]
    cdown=counts[,.(id1,id2,distance,count=contact.down)]
    data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
                 Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
                 Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
                 Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
                 Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
                 eC=par$eC, log_nu=par$log_nu, log_delta=par$log_delta,
                 beta_diag_centered=par$beta_diag_centered)
    as.data.table(optimizing(stanmodels$predict_all, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)$par)
  }
  return(pred)
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
run_split_parallel_initial_guess = function(counts, biases, bf_per_kb, dmin, dmax, bf_per_decade, verbose, iter) {
  #compute column sums
  cs1=counts[,.(R=sum(contact.close+contact.down),L=sum(contact.far+contact.up)),by=pos1][,.(pos=pos1,L,R)]
  cs2=counts[,.(R=sum(contact.far+contact.down),L=sum(contact.close+contact.up)),by=pos2][,.(pos=pos2,L,R)]
  pos=data.table(pos=positions, key="pos")
  sums=rbind(cs1,cs2)[,.(L=sum(L),R=sum(R)),keyby=pos][pos]
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  op=optimizing(stanmodels$guess, data=list(Krow=Krow, S=biases[,.N],
                                            cutsites=biases[,pos], rejoined=biases[,rejoined],
                                            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                            counts_sum_left=sums[,L], counts_sum_right=sums[,R]),
                as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
  #return initial guesses
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  return(list(eRJ=op$par$eRJ, eDE=op$par$eDE, eC=op$par$eC,
              alpha=op$par$alpha, beta_nu=op$par$beta_nu, beta_delta=op$par$beta_delta,
              lambda_nu=op$par$lambda_nu, lambda_delta=op$par$lambda_delta,
              log_nu=op$par$log_nu, log_delta=op$par$log_delta,
              beta_diag=rep(0,Kdiag-1), lambda_diag=1))
}

#' Genomic bias estimation part of parallel run
#' @keywords internal
#' 
run_split_parallel_biases = function(counts, biases, begin, end, dmin, dmax, bf_per_kb, bf_per_decade, verbose, iter, outprefix=NULL, circularize=-1) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases, begin1=begin, end1=end, fill.zeros=T, circularize=circularize)
  #run fit
  a=system.time(output <- capture.output(op <- csnorm_fit(stanmodels$fit, extracted$biases1, extracted$counts, dmin, dmax,
                                                          bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)))
  op$runtime=a[1]+a[4]
  op$output=output
  if (!is.null(outprefix)) {
    save(op, file=paste0(outprefix, "_biases_op_",begin,".RData"))
  }
  #report data
  center=(end-begin)/2
  ret=extracted$biases1[,.(id=id+extracted$beginrange1-1,pos,weight=(center-abs(pos-center))**2)]
  ret[,c("log_mean_DL", "log_mean_DR", "log_mean_RJ"):=list(op$par$log_mean_DL, op$par$log_mean_DR, op$par$log_mean_RJ)]
  ret[,square.begin:=begin]
  ret[,c("lambda_delta","lambda_nu","nsites","alpha"):=list(op$par$lambda_delta,op$par$lambda_nu,extracted$biases1[,.N],op$par$alpha)]
  if (!is.null(outprefix)) {
    save(ret, file=paste0(outprefix, "_biases_ret_",begin,".RData"))
  }
  return(list(ret=ret, out=tail(op$output,1), runtime=op$runtime))
}

#' Homogenize biases estimated on multiple fits
#' @keywords internal
#' 
run_split_parallel_biases_homogenize = function(dt, lambda_nu, lambda_delta, bf_per_kb=1, iter=10000, verbose=T, init=0, ...) {
  Krow=round(dt[,(max(pos)-min(pos))/1000*bf_per_kb])
  data = list( Krow=Krow, S=dt[,.N], cutsites=dt[,pos],
               log_mean_RJ=dt[,log_mean_RJ], log_mean_DL=dt[,log_mean_DL], log_mean_DR=dt[,log_mean_DR],
               lambda_nu=lambda_nu, lambda_delta=lambda_delta)
  optimizing(stanmodels$homogenize, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

#' Decay bias estimation part of parallel run
#' @keywords internal
#' 
run_split_parallel_counts = function(counts, biases.aug, begin1, end1, begin2, end2, dmin, dmax,
                                     bf_per_decade=5, verbose=T, iter=100000, outprefix=NULL, circularize=-1L) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases.aug, begin1=begin1, end1=end1,
                            begin2=begin2, end2=end2, fill.zeros=T, circularize=circularize)
  #run fit
  a=system.time(output <- capture.output(op <- csnorm_fit_extradiag(stanmodels$fit_extradiag, extracted$biases1, extracted$biases2, extracted$counts,
                                                                    dmin, dmax, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)))
  op$runtime=a[1]+a[4]
  op$output=output
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
  return(list(ret=ret, out=tail(op$output,1), runtime=op$runtime))
}

#' Merge together parallel runs and estimate a single diagonal decay
#' @keywords internal
#' 
run_split_parallel_counts_decay = function(dt, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  dist=dt[,(begin+end)/2]
  fij=dt[,(log_decay_close+log_decay_far+log_decay_up+log_decay_down)/4]
  ncounts=dt[,ncounts]
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=length(dist), dist=dist, fij=fij, ncounts=ncounts)
  optimizing(stanmodels$fit_decay, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

#' Get a proper count exposure for all estimated biases
#' @keywords internal
#' 
run_split_parallel_counts_eC = function(biases, counts, retlist, dmin, dmax, bf_per_decade=5, verbose=T, iter=1000) {
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
  optimizing(stanmodels$fit_eC, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
}


#' Parse intermediate outputs of killed run and try to compute the estimates nevertheless
#' @keywords internal
#' @export
#' 
run_split_parallel_recovery = function(cs, outprefix, square.size=100000, coverage=4, coverage.extradiag=1, bf_per_kb=1, bf_per_decade=5,
                                       distance_bins_per_decade=100, verbose = F, iter=100000, ncores=30, homogenize=F) {
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
  run_split_parallel(cs, NULL, square.size, coverage, coverage.extradiag, bf_per_kb, bf_per_decade,
                     distance_bins_per_decade, verbose, iter, ncores, homogenize, outprefix,
                     ops.bias, ops.count)
}

#' Plot bias run statistics using intermediate files
#' @keywords internal
#' @export
#' 
diagnose_biases = function(outprefix, coverage=4, square.size=150000) {
  diagsquares = run_split_parallel_squares(biases, square.size, coverage, diag.only=T)$diagsquares
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
  return(list(p0,p1,p2))
}


#' Plot bias run statistics using intermediate files
#' @keywords internal
#' @export
#' 
diagnose_counts = function(outprefix, coverage.extradiag=1, square.size=150000) {
  squares = run_split_parallel_squares(biases, square.size, coverage.extradiag, diag.only=F)$squares
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
  return(list(p0,p1,p2))
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
#' @return A list containing:
#' \enumerate{
#' \item par: the optimized parameters
#' \item out.bias: the last line of the stan output for the bias estimation
#' \item runtime.bias: the runtimes of these optimizations
#' \item out.count: the last line of the stan output for the decay estimation
#' \item runtime.count: the runtimes of these optimizations
#' }
#' @export
#' 
#' @examples
run_split_parallel = function(cs, design=NULL, square.size=100000, coverage=4, coverage.extradiag=1, bf_per_kb=1, bf_per_decade=5,
                              distance_bins_per_decade=100, verbose = F, iter=100000, ncores=30, homogenize=F, outprefix=NULL,
                              ops.bias=NULL, ops.count=NULL) {
  stopifnot(!(cs@settings$circularize==-1 && cs@counts[,max(distance)]<cs@biases[,max(pos)-min(pos)]))
  cs@settings = c(cs@settings, list(square.size=square.size, coverage=coverage, coverage.extradiag=coverage.extradiag,
                                    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, distance_bins_per_decade=distance_bins_per_decade,
                                    iter=iter, homogenize=homogenize, outprefix=outprefix))
  stopifnot(is.null(design))
  ### build squares
  message("*** build squares")
  retval.squares = run_split_parallel_squares(cs@biases, square.size, coverage, diag.only=T)
  diagsquares=retval.squares$diagsquares
  true.square.size=retval.squares$true.square.size
  retval.squares = run_split_parallel_squares(cs@biases, square.size, coverage.extradiag, diag.only=F)
  stopifnot(true.square.size==retval.squares$true.square.size)
  squares=retval.squares$squares
  message("Run in parallel")
  message("   square size: ",true.square.size)
  message("   diagonal coverage: ",coverage)
  message("   extradiagonal coverage: ",coverage.extradiag)
  message("   ",squares[,.N], " squares")
  message("   ",diagsquares[,.N], " on the diagonal")
  ### fit genomic biases using squares close to diagonal
  message("*** fit genomic biases")
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
    list(par=rbindlist(a[,1]), out=as.character(a[,2]), runtime=as.numeric(a[,3]))
  }
  if (is.null(ops.bias)) {
    ops.bias = foreach (begin=diagsquares[,begin1], end=diagsquares[,end1], .packages=c("data.table","rstan")) %dopar% 
      run_split_parallel_biases(cs@counts, cs@biases, begin, end, dmin, dmax, bf_per_kb = bf_per_kb,
                                bf_per_decade = bf_per_decade, verbose = verbose, iter = iter, outprefix=outprefix,
                                circularize=cs@settings$circularize)
    ops.bias = output.binder(ops.bias)
    if (!is.null(outprefix)) {save(ops.bias, file=paste0(outprefix, "_ops_bias.RData"))}
  }
  ### reconstruct bias estimates: eRJ eDE log_nu and log_delta
  message("*** reconstruct genomic biases")
  info=ops.bias$par[,.SD[1],by=square.begin,.SDcols=c("lambda_nu","lambda_delta","nsites","alpha")]
  #ggplot(info)+geom_line(aes(square.begin,lambda_nu))+scale_y_log10()+geom_hline(yintercept=info[,exp(median(log(lambda_nu)))])
  means=ops.bias$par[,.(log_mean_RJ=weighted.mean(log_mean_RJ, weight),
                        log_mean_DL=weighted.mean(log_mean_DL, weight),
                        log_mean_DR=weighted.mean(log_mean_DR, weight),
                        ncounts=.N), keyby=c("id","pos")]
  if (homogenize==T) {
    message("*** homogenize genomic biases")
    op=run_split_parallel_biases_homogenize(means, lambda_nu=info[,exp(mean(log(lambda_nu)))], lambda_delta=info[,exp(mean(log(lambda_delta)))],
                                            bf_per_kb=bf_per_kb, iter=iter, verbose=verbose)
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
  #ggplot()+geom_line(data=data.table(x=biases.aug[,pos],y=retlist$log_nu),aes(x,y))
  #
  ### fit remaining data
  message("*** fit diagonal decay")
  registerDoParallel(cores=ncores)
  if (is.null(ops.count)) {
    ops.count = foreach (begin1=squares[,begin1], end1=squares[,end1], begin2=squares[,begin2], end2=squares[,end2],
                         .packages=c("data.table","rstan")) %dopar% 
      run_split_parallel_counts(cs@counts, biases.aug, begin1, end1, begin2, end2, dmin, dmax,
                                bf_per_decade = bf_per_decade, verbose = verbose, iter = iter, outprefix=outprefix,
                                circularize=cs@settings$circularize)
    ops.count=output.binder(ops.count)
    if (!is.null(outprefix)) {save(ops.count, file=paste0(outprefix, "_ops_count.RData"))}
  }
  ### reconstruct count estimates: log_decay and beta_diag_centered
  message("*** reconstruct diagonal decay")
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
  op=run_split_parallel_counts_decay(test, counts, dmin, dmax, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)
  test[,log_decay:=op$par$log_decay+op$par$eC]
  retlist$decay=test[,.(dist=(begin+end)/2,decay=exp(log_decay))]
  retlist$beta_diag_centered=op$par$beta_diag_centered
  ### reconstruct count estimates: eC
  message("*** reconstruct count exposure")
  op=run_split_parallel_counts_eC(cs@biases, cs@counts[sample(.N,min(.N,100000))], retlist, dmin, dmax, bf_per_decade=bf_per_decade, verbose=verbose)
  retlist$eC=op$par$eC
  retlist$alpha=op$par$alpha
  #ggplot(test)+geom_point(aes(bdist, log_decay))
  #ggplot(melt(test, id.vars=c("bdist","begin","end")))+geom_point(aes(bdist,value,colour=variable))
  message("*** predict all means")
  pred=csnorm_predict_all(cs, verbose=verbose, ncores=ncores)
  pred #for some reason it's needed
  cs@pred=pred
  message("*** done")
  cs@par=retlist
  cs@diagnostics=list(out.bias=ops.bias$out, out.count=ops.count$out, runtime.count=ops.count$runtime, runtime.bias=ops.bias$runtime)
  return(cs)
}
