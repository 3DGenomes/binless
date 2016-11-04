#' @include csnorm.R
NULL

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#' 
csnorm_gauss_genomic_bam = function(cs, verbose=T, init.mean="mean", nthreads=1) {
  if (init.mean=="mean") {
    a = csnorm:::csnorm_gauss_genomic_muhat_mean(cs)
  } else {
    a = csnorm_gauss_genomic_muhat_data(cs)
  }
  bts=a$bts
  cts=a$cts
  #add factors for design matrix
  data=rbind(bts,cts)
  data[,weight:=1/std^2]
  data[cat=="contact L",  c("iota.coef","rho.coef","is.dangling","is.rejoined"):=list(1,0,0,0)]
  data[cat=="contact R",  c("iota.coef","rho.coef","is.dangling","is.rejoined"):=list(0,1,0,0)]
  data[cat=="dangling L", c("iota.coef","rho.coef","is.dangling","is.rejoined"):=list(1,0,1,0)]
  data[cat=="dangling R", c("iota.coef","rho.coef","is.dangling","is.rejoined"):=list(0,1,1,0)]
  data[cat=="rejoined",   c("iota.coef","rho.coef","is.dangling","is.rejoined"):=list(1/2,1/2,0,1)]
  Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
  #optimize each set of biases separately
  cat("optimize with BAM\n") #necessary so logging goes smoothly
  sdata = foreach(gen=cs@design[,uniqueN(genomic)],
                 .combine=function(x,y){list(eta=rbind(x$eta,y$eta),value=x$value+y$value,sp=c(x$sp,y$sp))}) %do% {
    dsets=cs@design[genomic==gen,name]
    sdata=data[name %in% dsets]
    #formula differs if there is more than one dataset
    if (length(dsets)>1) {
      sdata[,dset:=as.integer(factor(name))-1]
      fit=bam(formula = etahat ~ dset+is.dangling+is.rejoined-1
                         +s(pos, by=iota.coef,bs="ps", m=2, k=Krow) +s(pos, by=rho.coef,bs="ps", m=2, k=Krow),
              data=sdata, family=gaussian(), weight=sdata[,weight], scale=1, discrete=T, samfrac=0.1, nthreads=nthreads)
    } else {
      fit=bam(formula = etahat ~ is.dangling+is.rejoined-1
              +s(pos, by=iota.coef,bs="ps", m=2, k=Krow) +s(pos, by=rho.coef,bs="ps", m=2, k=Krow),
              data=sdata, family=gaussian(), weight=sdata[,weight], scale=1, discrete=T, samfrac=0.1, nthreads=nthreads)
    }
    list(eta=fit$fitted.values,value=fit$gcv.ubre,sp=list(fit$sp))
  }
  #return values as if optimized by stan
  data[,eta:=sdata$eta]
  setkeyv(data,key(cs@biases))
  eC=data[cat%in%c("contact L","contact R"), .(eC=mean(eta)),by=name]
  eDE=data[cat%in%c("dangling L","dangling R"), .(eDE=mean(eta)),by=name]
  eRJ=data[cat=="rejoined", .(eRJ=mean(eta)),by=name]
  log_iota=data[cat=="contact L",.(id,log_iota=eta-mean(eta)),by=name]
  log_rho=data[cat=="contact R",.(id,log_rho=eta-mean(eta)),by=name]
  data=data[merge(log_iota,log_rho,by=key(data))]
  lfac = 30*cs@settings$bf_per_kb
  op=list(value=sdata$value, par=list(eC=eC[,as.array(eC)], eDE=eDE[,as.array(eDE)], eRJ=eRJ[,as.array(eRJ)],
                             log_iota=log_iota[,log_iota], log_rho=log_rho[,log_rho],
                             lambda_iota=sqrt(as.array(sapply(sdata$sp,function(x){x[[1]]})))/lfac,
                             lambda_rho=sqrt(as.array(sapply(sdata$sp,function(x){x[[2]]})))/lfac))
  op$par$biases = data[,.(cat, name, id, pos, etahat, std, eta)]
  setkey(op$par$biases, id, name, cat)
  #update par slot
  op$par$value=op$value
  cs@par=modifyList(cs@par, op$par)
  return(cs)
}

#' Cut-site normalization (deprecated fast approximation, using mgcv BAM)
#'
#' Note that the genomic length scales are optimized twice during the
#' optimization
#' 
#' @inheritParams run_gauss
#' @export
#' 
#' @examples
run_gauss_bam = function(cs, init=NULL, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10,
                           ngibbs = 3, iter=10000, fit.decay=T, fit.genomic=T, fit.disp=T,
                           verbose=T, ncounts=100000, init_alpha=1e-7, ncores=1) {
  #clean object if dirty
  cs@par=list() #in case we have a weird object
  cs@binned=list()
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
                                    iter=iter, ngibbs=ngibbs, init_alpha=init_alpha))
  #fill counts matrix and take subset
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize, dmin=cs@settings$dmin)
  setkey(cs@biases, id, name)
  setkey(cs@counts, id1, id2, name)
  ncounts_per_dset=as.integer(ncounts/cs@design[,.N])
  subcounts = cs@counts[,.SD[1:min(.N,ncounts_per_dset)],by=name]
  setkeyv(subcounts,key(cs@counts))
  if (subcounts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2)
    stop("dataset too sparse, please increase ncounts")
  subcounts.weight=merge(cs@counts[,.(nc=.N),keyby=name],subcounts[,.(ns=.N),keyby=name])[,.(name,wt=nc/ns)]
  #initial guess
  if (is.null(init)) {
    if (verbose==T) cat("No initial guess provided\n")
    cs@diagnostics=list()
    laststep=0
    init.mean="data"
    cs=fill_parameters(cs, dispersion=1, fit.decay=fit.decay, fit.genomic=fit.genomic, fit.disp=fit.disp) #init with dispersion=1
  } else {
    if (verbose==T) cat("Using provided initial guess\n")
    if (is.data.table(cs@diagnostics$params)) laststep = cs@diagnostics$params[,max(step)]
    init$beta_diag = guarantee_beta_diag_increasing(init$beta_diag)
    init.mean="mean"
    cs@par=init
  }
  #gibbs sampling
  for (i in (laststep + 1:ngibbs)) {
    #fit diagonal decay given iota and rho
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay\n")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_decay(cs, init.mean=init.mean, init_alpha=init_alpha)))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="decay", out=output, runtime=a[1]+a[4])
    }
    #fit iota and rho given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic\n")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_genomic_bam(cs, init.mean=init.mean, nthreads=ncores)))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="bias", out=output, runtime=a[1]+a[4])
    }
    init.mean="mean"
    if (fit.disp==T) {
      #fit exposures and dispersion
      if (verbose==T) cat("Gibbs",i,": Remaining parameters\n")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_dispersion(cs, counts=subcounts, weight=subcounts.weight,
                                                                                    init_alpha=init_alpha)))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="disp", out=output, runtime=a[1]+a[4])
      if (verbose==T) cat("Gibbs",i,": log-likelihood = ",cs@par$value,"\n")
    }
  }
  if (verbose==T) cat("Done\n")
  cs@par$init=init
  cs@diagnostics$plot=ggplot(cs@diagnostics$params[,.(step,leg,value,out.last)])+
    geom_line(aes(step,value))+geom_point(aes(step,value,colour=out.last))+facet_wrap(~leg, scales = "free")
  return(cs)
}

