#' Bin a counts data.table into a matrix of a given resolution
#'
#' @param counts data.table as returned by \code{\link{prepare_for_sparse_cs_norm}}
#' @param resolution positive integer.
#' @param b1, b2, e1, e2 Begins and ends of the portion of the data to bin. If NULL replace by min/max value. 
#'
#' @return a data.table representing the binned data
#' @export
#'
#' @examples
bin_counts = function(counts, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL) {
  if (is.null(b1)) b1=counts[,min(pos1)]-1
  if (is.null(b2)) b2=counts[,min(pos2)]-1
  if (is.null(e1)) e1=counts[,max(pos1)]-1
  if (is.null(e2)) e2=counts[,max(pos2)]-1
  bins1=seq(b1,counts[,max(pos1)]+resolution,resolution)
  bins2=seq(b2,counts[,max(pos2)]+resolution,resolution)
  mcounts=melt(counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0]
  #
  sub = mcounts[,.(pos1,pos2,bin1=cut2(pos1, bins1, oneval=F, onlycuts=T, digits=10),
                   bin2=cut2(pos2, bins2, oneval=F, onlycuts=T, digits=10), category, count)
                ][,.(N=sum(count)),by=c("bin1","bin2")]
  #
  sub[,begin1:=do.call(as.integer, tstrsplit(as.character(bin1), "[[,]")[2])]
  sub[,end1:=do.call(as.integer, tstrsplit(as.character(bin1), "[],)]")[2])]
  sub[,begin2:=do.call(as.integer, tstrsplit(as.character(bin2), "[[,]")[2])]
  sub[,end2:=do.call(as.integer, tstrsplit(as.character(bin2), "[],)]")[2])]
  #
  return(sub)
}

#' Apply the ICE algorithm to a binned matrix
#'
#' @param bdata binned data.table as returned by \code{\link{bin_counts}}
#' @param niterations positive integer. Number of iterations to perform
#'
#' @return the ICEd matrix
#' @export
#'
#' @examples
iterative_normalization = function(bdata, niterations=100) {
  binned = bdata[bin1<bin2,.(bin1,bin2,N=observed)]
  binned = rbind(binned, binned[,.(bin1=bin2,bin2=bin1,N)])
  binned[,N.weighted:=N]
  #iterate
  for (i in 1:niterations) {
    binned[,b1:=sum(N.weighted),by=bin1]
    binned[,b2:=sum(N.weighted),by=bin2]
    binned[,b1:=b1/mean(b1)]
    binned[,b2:=b2/mean(b2)]
    binned[,N.weighted:=N.weighted/b1/b2]
  }
  binned[,c("b1","b2","N"):=list(NULL,NULL,NULL)]
  setnames(binned,"N.weighted","N")
  binned[,begin1:=do.call(as.integer, tstrsplit(as.character(bin1), "[[,]")[2])]
  binned[,end1:=do.call(as.integer, tstrsplit(as.character(bin1), "[],)]")[2])]
  binned[,begin2:=do.call(as.integer, tstrsplit(as.character(bin2), "[[,]")[2])]
  binned[,end2:=do.call(as.integer, tstrsplit(as.character(bin2), "[],)]")[2])]
  return(binned[begin1<begin2])
}


#' Predict expected values for each count given optimized model parameters
#'
#' @param model the stan model
#' @param biases, counts as returned by \code{\link{prepare_for_sparse_cs_norm}} or a subset of it
#' @param opt the optimized parameters, as returned by \code{\link{run_split_parallel}} 
#' @param bf_per_decade same as that used in \code{\link{run_split_parallel}}
#' @param verbose 
#'
#' @return a stan optimization output with the predictions.
#' @keywords internal
#' @export
#'
#' @examples
csnorm_predict_all = function(model, biases, counts, opt, bf_per_decade=5, verbose=T) {
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round((log10(opt$dmax)-log10(opt$dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=opt$dmin, dmax=opt$dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag_centered=opt$par$beta_diag_centered)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
}

#' Bin observed and expected counts at a given resolution
#' @keywords internal
#' @export
#' 
csnorm_predict_binned = function(model, biases, counts, opt, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL,
                               bf_per_decade=5, verbose=F, circularize=-1L) {
  csub=copy(counts) #need to implement taking only needed part of matrix, and reporting log_nu and log_delta appropriately
  bsub=copy(biases)
  Kdiag=round((log10(opt$dmax)-log10(opt$dmin))*bf_per_decade)
  npoints=100*Kdiag #evaluate spline with 100 equidistant points per basis function
  #bin existing counts and biases
  if (is.null(b1)) b1=bsub[,min(pos)]-1
  if (is.null(b2)) b2=bsub[,min(pos)]-1
  if (is.null(e1)) e1=bsub[,max(pos)]+1
  if (is.null(e2)) e2=bsub[,max(pos)]+1
  bins1=seq(b1,e1+resolution,resolution)
  bins2=seq(b2,e2+resolution,resolution)
  csub[,c("bin1","bin2"):=list(cut(pos1, bins1, ordered_result=T, right=F, include.lowest=T,dig.lab=5),
                               cut(pos2, bins2, ordered_result=T, right=F, include.lowest=T,dig.lab=5))]
  csub[,log_mean_cup:=opt$pred$log_mean_cup]
  csub[,log_mean_cdown:=opt$pred$log_mean_cdown]
  csub[,log_mean_cclose:=opt$pred$log_mean_cclose]
  csub[,log_mean_cfar:=opt$pred$log_mean_cfar]
  bsub[,log_nu:=opt$par$log_nu]
  bsub[,log_delta:=opt$par$log_delta]
  bsub[,c("bin1","bin2"):=list(cut(pos, bins1, ordered_result=T, right=F, include.lowest=T,dig.lab=5),
                               cut(pos, bins2, ordered_result=T, right=F, include.lowest=T,dig.lab=5))]
  bsub[,c("ibin1","ibin2"):=list(as.integer(bin1),as.integer(bin2))]
  #run computation of mean
  cclose=csub[,.(id1,id2,bin1,bin2,distance,count=contact.close,logmean=log_mean_cclose)]
  cfar=csub[,.(id1,id2,bin1,bin2,distance,count=contact.far,logmean=log_mean_cfar)]
  cup=csub[,.(id1,id2,bin1,bin2,distance,count=contact.up,logmean=log_mean_cup)]
  cdown=csub[,.(id1,id2,bin1,bin2,distance,count=contact.down,logmean=log_mean_cdown)]
  csub=rbind(cclose,cfar,cup,cdown)[!(is.na(bin1)|is.na(bin2)|count==0)]
  data = list( Kdiag=Kdiag, npoints=npoints, circularize=circularize,
               S1=bsub[!is.na(bin1),.N], S2=bsub[!is.na(bin2),.N], 
               cutsites1=bsub[!is.na(bin1),pos], cutsites2=bsub[!is.na(bin2),pos],
               dmin=opt$dmin, dmax=opt$dmax,
               N=csub[,.N],   counts=csub[,count], cdist=csub[,distance], cmean=csub[,exp(logmean)],
               eC=opt$par$eC,
               log_nu1=bsub[!is.na(bin1),log_nu], log_nu2=bsub[!is.na(bin2),log_nu],
               log_delta1=bsub[!is.na(bin1),log_delta], log_delta2=bsub[!is.na(bin2),log_delta],
               beta_diag_centered=opt$par$beta_diag_centered,
               B1=csub[,nlevels(bin1)], B2=csub[,nlevels(bin2)],
               cbins1=csub[,as.integer(bin1)], cbins2=csub[,as.integer(bin2)],
               bbins1=bsub[!is.na(bin1),ibin1], bbins2=bsub[!is.na(bin2),ibin2])
  binned=optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
  #format output
  so = data.table(melt(binned$par$observed))
  setnames(so,"value","observed")
  so[,expected:=melt(binned$par$expected)$value]
  so[,ncounts:=melt(binned$par$ncounts)$value]
  so[,normalized:=melt(binned$par$normalized)$value]
  so=so[ncounts>0]
  setkey(so,Var1)
  so=bsub[,.(bin1=bin1[1]),keyby=ibin1][so]
  setkey(so,Var2)
  so=bsub[,.(bin2=bin2[1]),keyby=ibin2][so]
  so[,c("ibin1","ibin2"):=list(NULL,NULL)]
  so[,lFC:=log2(observed/expected)]
  return(list(distance=exp(binned$par$log_dist), decay=exp(binned$par$log_decay), mat=so))
}

#' Predict dispersions for a binned matrix
#' @keywords internal
#' @export
#' 
get_dispersions = function(model, binned, iter=10000, verbose=T) {
  data=list(B=binned[,.N],observed=binned[,observed],expected=binned[,expected],ncounts=binned[,ncounts])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
}

#' compute \f$p(\Gamma_2>\Gamma_1) = \int_{0}^{+\infty} dx p_{\Gamma_2}(x) \int_{0}^{x} dy p_{\Gamma_1}(y)\f$
#' @keywords internal
#' 
compute_gamma_overlap = function(alpha1,beta1,alpha2,beta2, bounds=5, ncores=1) {
  registerDoParallel(cores=ncores)
  foreach (a1=alpha1, b1=beta1, a2=alpha2, b2=beta2, .packages="stats", .combine=c) %dopar% {
    #infinite integral does not work very well so truncate outer integral
    mu1=a1/b1
    mu2=a2/b2
    sd1=sqrt(a1)/b1
    sd2=sqrt(a2)/b2
    xmin = max(0,min(mu1-bounds*sd1,mu2-bounds*sd2))
    xmax = max(mu1+bounds*sd1,mu2+bounds*sd2)
    a=integrate(function(x){exp(dgamma(x,a2,rate=b2,log=T)+pgamma(x,a1,rate=b1,log.p=T))},xmin,xmax, stop.on.error=F)
    if (a$abs.error<=0 | a$message != "OK") {NA} else {a$value}
  }
}

#' compute \f$p(\mathcal{N}_2>\mathcal{N}_1) = \int_{-infty}^{+infty} dx p_{\mathcal{N}_2}(x) \int_{-infty}^{x} dy p_{\mathcal{N}_1}(y)\f$
#' @keywords internal
#' 
compute_normal_overlap = function(mu1,sd1,mu2,sd2, bounds=5, ncores=1) {
  registerDoParallel(cores=ncores)
  foreach (m1=mu1, s1=sd1, m2=mu2, s2=sd2, .packages="stats", .combine=c) %dopar% {
    #infinite integral does not work very well so truncate outer integral
    xmin = min(m1-bounds*s1,m2-bounds*s2)
    xmax = max(m1+bounds*s1,m2+bounds*s2)
    a=integrate(function(x){exp(dnorm(x,mean=(m2-m1)/s1,sd=s2/s1,log=T)+pnorm(x,mean=0,sd=1,log.p=T))},(xmin-m1)/s1,(xmax-m1)/s1)
    if (a$abs.error<=0 | a$message != "OK") {NA} else {a$value}
  }
}

#' Detect significant interactions wrt expected
#'
#' @param binned as returned by \code{\link{csnorm_predict_binned}}
#' @param dispersions as returned by \code{\link{get_dispersions}}
#' @param threshold significance threshold, between 0 and 1
#' @param ncores number of cores used for parallelization
#' @param normal.approx use normal approximation if dispersion is larger than this
#'
#' @return the binned matrix with additional information relating to these significant interactions
#' @keywords internal
#' @export
#'
#' @examples
detect_interactions = function(binned, dispersions, threshold=0.95, ncores=1, normal.approx=100){
  #report gamma parameters
  mat=copy(binned)
  mat[,c("alpha1","beta1"):=list(dispersions,dispersions/expected)]
  mat[,c("alpha2","beta2"):=list(alpha1+observed,beta1+1)]
  mat[alpha1<normal.approx,c("prob.observed.gt.expected","detection.type"):=list(compute_gamma_overlap(alpha1,beta1,alpha2,beta2,ncores=ncores),"gamma")]
  mat[,c("mean1","sd1"):=list(expected,expected/sqrt(dispersions))]
  mat[,c("mean2","sd2"):=list(alpha2/beta2, sqrt(alpha2)/beta2)]
  mat[alpha1>=normal.approx,c("prob.observed.gt.expected","detection.type"):=list(compute_normal_overlap(mean1,sd1, mean2, sd2, ncores=ncores),"normal")]
  mat[,prob.observed.gt.expected:=as.numeric(prob.observed.gt.expected)]
  mat[is.na(prob.observed.gt.expected)&observed>=expected,c("prob.observed.gt.expected","detection.type"):=list(ppois(observed,expected,lower.tail=F),"poisson")]
  mat[is.na(prob.observed.gt.expected)&observed<expected,c("prob.observed.gt.expected","detection.type"):=list(ppois(observed,expected,lower.tail=T),"poisson")]
  mat[,is.interaction:=prob.observed.gt.expected>threshold | 1-prob.observed.gt.expected>threshold]
  #write begins/ends
  bin1.begin=mat[,bin1]
  bin1.end=mat[,bin1]
  bin2.begin=mat[,bin2]
  bin2.end=mat[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[[,]")[2][[1]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[[,)]")[2][[1]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[[,]")[2][[1]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[[,)]")[2][[1]]
  mat[,begin1:=as.integer(as.character(bin1.begin))]
  mat[,end1:=as.integer(as.character(bin1.end))]
  mat[,begin2:=as.integer(as.character(bin2.begin))]
  mat[,end2:=as.integer(as.character(bin2.end))]
  return(mat)
}

#' estimates the values of the count or dispersion required to cross a given threshold
#' 
#' For illustration purposes. Has a border effect at high counts that can be safely ignored.
#'
#' @param observed, expected, dispersion, threshold floats used to generate the plots 
#' @param counts.range scan counts in that range
#' @param disp.range scan dispersion in that range
#' @param compute whether to compute the values that meet the threshold or not. This step can fail if the ranges are set badly.
#'
#' @return a list of three plots
#' @export
#'
#' @examples
thresholds_estimator = function(observed, expected, dispersion, threshold=0.95, counts.range=c(expected/100,expected*100), disp.range=c(0,100*dispersion), compute=T) {
  #plot gamma distributions
  a1=dispersion
  b1=dispersion/expected
  a2=a1+observed
  b2=b1+1
  mu1=a1/b1
  mu2=a2/b2
  sd1=sqrt(a1)/b1
  sd2=sqrt(a2)/b2
  message("prior:     mean ",mu1, " sd ",sd1)
  message("posterior: mean ",mu2, " sd ",sd2)
  xmin = min(c(0,counts.range))
  xmax = max(counts.range)
  message("evaluation bounds: [",xmin,",",xmax,"]")
  p1=ggplot(data.table(x=c(xmin,xmax)),aes(x))+
    stat_function(fun=function(x){dgamma(x,a1,rate=b1)},aes(colour="prior"))+
    stat_function(fun=function(x){dgamma(x,a2,rate=b2)},aes(colour="posterior"))+
    scale_colour_manual("Gamma", values = c("blue", "red"))
  #helper function  
  gammaQuery=function(mu,theta,count){
    alpha1=theta
    alpha2=theta+count
    beta1=theta/mu
    beta2=beta1+1
    return(function(x){exp(dgamma(x,alpha2,rate=beta2,log=T)+pgamma(x,alpha1,rate=beta1,log.p=T))})
  }
  #query counts
  igammaQuery_c=function(x){integrate(gammaQuery(expected,dispersion,count=x),xmin,xmax)$value-threshold}
  p2=ggplot(data.table(x=c(xmin,xmax),y=c(0,1)),aes(x))+
    stat_function(fun=Vectorize(function(x){igammaQuery_c(x)+threshold}))+
    geom_hline(aes(yintercept=threshold))+ geom_hline(aes(yintercept=1-threshold))+
    labs(x="observed count", y="P(observed>expected)")
  if (compute){
    clo = uniroot(function(x){igammaQuery_c(x)+2*threshold-1},c(xmin,expected))$root
    chi = uniroot(igammaQuery_c,c(expected,xmax))$root
    message("At a threshold of ",threshold," significant counts are lower than ",clo," and higher than ",chi)
    p2=p2+geom_vline(aes(xintercept=clo))+ geom_vline(aes(xintercept=chi))
  }
  #query dispersion
  igammaQuery_theta=function(x){integrate(gammaQuery(expected,x,observed),xmin,xmax)$value-threshold}
  p3=ggplot(data.table(x=disp.range,y=c(0,1)),aes(x))+
    stat_function(fun=Vectorize(function(x){igammaQuery_theta(x)+threshold}))+
    geom_hline(aes(yintercept=threshold))+geom_hline(aes(yintercept=1-threshold))+
    labs(x="dispersion", y="P(observed>expected)")+scale_x_log10()
  if (compute) {
    dstar = uniroot(igammaQuery_theta,disp.range)$root
    p3=p3+geom_vline(aes(xintercept=dstar))
    if (observed > expected) {
      message("At a threshold of ",threshold," the dispersion must be larger than ",dstar," to reach significance")
    } else {
      message("At a threshold of ",threshold," the dispersion must be smaller than ",dstar," to reach significance")
    }
  }
  return(list(p1,p2,p3))
}

#' Wrapper for the postprocessing steps
#'
#' @param biases 
#' @param counts 
#' @param op 
#' @param resolution 
#' @param ncores 
#' @param predict.all.means 
#' @param circularize 
#'
#' @return
#' @export
#'
#' @examples
postprocess = function(biases, counts, op, resolution=10000, ncores=30, predict.all.means=T, circularize=-1L) {
  smpred = stan_model(file = "predict_all.stan")
  smbin = stan_model("predict_binned.stan")
  smdisp = stan_model("dispersions.stan")
  ### run remaining steps
  if (predict.all.means==T) {
    message("*** predict all means")
    op$pred=csnorm_predict_all(smpred, biases, counts, op, verbose=T)$par
  }
  message("*** buid binned matrices")
  op$binned=csnorm_predict_binned(smbin, biases, counts, op, resolution=resolution, circularize=circularize)
  message("*** estimate dispersions")
  op$disp=get_dispersions(smdisp, op$binned$mat)$par
  message("*** detect interactions")
  op$mat=detect_interactions(op$binned$mat, op$disp$dispersion, ncores=ncores) #interaction detection using binned dispersion estimates
  return(op)
}

