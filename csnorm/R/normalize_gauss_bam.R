#' @include csnorm.R
NULL

#' Single-cpu simplified initial guess
#' @keywords internal
#' 
csnorm_gauss_guess_bam = function(biases, counts, design, dmin, dmax, bf_per_kb=1, bf_per_decade=20,
                              dispersion=10, nthreads=1, verbose=T) {
  op = csnorm_gauss_genomic_bam(biases, counts, design, init=dispersion,
                                bf_per_kb=bf_per_kb, verbose=verbose, nthreads=nthreads)
  #add diagonal decay inits
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  op$par=c(list(beta_diag=beta_diag, lambda_diag=array(1,dim=Decays), log_decay=rep(0,counts[,.N]),
                alpha=dispersion),
           op$par[c("eC","eRJ","eDE","log_iota","log_rho","biases","lambda_iota","lambda_rho")])
  op
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#' 
csnorm_gauss_genomic_bam = function(biases, counts, design, init, bf_per_kb=1, verbose=T, ...) {
  if (length(init)>1) {
    a = csnorm_gauss_genomic_muhat_mean(biases, counts, design, init)
  } else {
    a = csnorm_gauss_genomic_muhat_data(biases, counts, init[[1]])
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
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  #optimize each set of biases separately
  cat("optimize with BAM\n") #necessary so logging goes smoothly
  sdata = foreach(gen=design[,uniqueN(genomic)],
                 .combine=function(x,y){list(eta=rbind(x$eta,y$eta),value=x$value+y$value,sp=c(x$sp,y$sp))}) %do% {
    dsets=design[genomic==gen,name]
    sdata=data[name %in% dsets]
    #formula differs if there is more than one dataset
    if (length(dsets)>1) {
      sdata[,dset:=as.integer(factor(name))-1]
      fit=bam(formula = etahat ~ dset+is.dangling+is.rejoined-1
                         +s(pos, by=iota.coef,bs="ps", m=2, k=Krow) +s(pos, by=rho.coef,bs="ps", m=2, k=Krow),
              data=sdata, family=gaussian(), weight=sdata[,weight], scale=1, discrete=T, samfrac=0.1, ...)
    } else {
      fit=bam(formula = etahat ~ is.dangling+is.rejoined-1
              +s(pos, by=iota.coef,bs="ps", m=2, k=Krow) +s(pos, by=rho.coef,bs="ps", m=2, k=Krow),
              data=sdata, family=gaussian(), weight=sdata[,weight], scale=1, discrete=T, samfrac=0.1, ...)
    }
    list(eta=fit$fitted.values,value=fit$gcv.ubre,sp=list(fit$sp))
  }
  #return values as if optimized by stan
  data[,eta:=sdata$eta]
  setkeyv(data,key(biases))
  eC=data[cat%in%c("contact L","contact R"), .(eC=mean(eta)),by=name]
  eDE=data[cat%in%c("dangling L","dangling R"), .(eDE=mean(eta)),by=name]
  eRJ=data[cat=="rejoined", .(eRJ=mean(eta)),by=name]
  log_iota=data[cat=="contact L",.(id,log_iota=eta-mean(eta)),by=name]
  log_rho=data[cat=="contact R",.(id,log_rho=eta-mean(eta)),by=name]
  data=data[merge(log_iota,log_rho,by=key(data))]
  lfac = 30*bf_per_kb
  op=list(value=sdata$value, par=list(eC=eC[,as.array(eC)], eDE=eDE[,as.array(eDE)], eRJ=eRJ[,as.array(eRJ)],
                             log_iota=log_iota[,log_iota], log_rho=log_rho[,log_rho],
                             lambda_iota=sqrt(as.array(sapply(sdata$sp,function(x){x[[1]]})))/lfac,
                             lambda_rho=sqrt(as.array(sapply(sdata$sp,function(x){x[[2]]})))/lfac))
  op$par$biases = data[,.(cat, name, id, pos, etahat, std, eta)]
  setkey(op$par$biases, cat, name, id)
  op
}

#' Cut-site normalization (deprecated fast approximation, using mgcv BAM)
#' 
#' @inheritParams run_gauss
#' 
#' @examples
run_gauss_bam = function(cs, init=NULL, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10,
                           ngibbs = 3, iter=100000, fit.decay=T, fit.genomic=T, fit.disp=T,
                           verbose=T, ncounts=100000, init_alpha=1e-5, ncores=1) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
                                    iter=iter, ngibbs=ngibbs))
  #fill counts matrix and take subset
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize, dmin=cs@settings$dmin)
  setkey(cs@biases, name, id)
  setkey(cs@counts, name, id1, id2)
  ncounts_per_dset=as.integer(ncounts/cs@design[,.N])
  subcounts = cs@counts[,.SD[1:min(.N,ncounts_per_dset)],by=name]
  if (subcounts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2)
    stop("dataset too sparse, please increase ncounts")
  subcounts.weight=merge(cs@counts[,.(nc=.N),keyby=name],subcounts[,.(ns=.N),keyby=name])[,.(name,wt=nc/ns)]
  dmin=cs@settings$dmin
  dmax=cs@settings$dmax
  #initial guess
  if (is.null(init)) {
    if (verbose==T) cat("No initial guess provided\n")
    cs@par=list() #in case we have a weird object
    cs@diagnostics=list()
    cs@binned=list()
    laststep=0
    init.a = system.time(init.output <- capture.output(op <- csnorm:::csnorm_gauss_guess_bam(
      biases = cs@biases, counts = cs@counts, design = cs@design, dmin=dmin, dmax=dmax,
      bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, dispersion=1, nthreads=ncores)))
    init.output = "Init with data"
    cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=0, leg="bias", out=init.output,
                                                 runtime=init.a[1]+init.a[4], op=op)
  } else {
    if (verbose==T) cat("Using provided initial guess\n")
    init.output = "Init provided"
    op=list(value=-1, par=init)
    cs@binned=list()
    if (is.data.table(cs@diagnostics$params)) laststep = cs@diagnostics$params[,max(step)]
  }
  #make sure beta_diag is strictly increasing
  op$par$beta_diag = guarantee_beta_diag_increasing(op$par$beta_diag)
  #gibbs sampling
  for (i in (laststep + 1:ngibbs)) {
    #fit diagonal decay given iota and rho
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay\n")
      a=system.time(output <- capture.output(op.diag <- csnorm:::csnorm_gauss_decay(
        biases = cs@biases, counts = cs@counts, design=cs@design,
        init=op$par, dmin = dmin, dmax = dmax,
        bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter, init_alpha=init_alpha)))
      op=list(value=op.diag$value, par=c(op.diag$par[c("eC","beta_diag","beta_diag_centered",
                                                       "log_decay","decay", "lambda_diag")],
                                         op$par[c("eRJ","eDE", "alpha",
                                                  "lambda_iota","lambda_rho","log_iota","log_rho","biases")]))
      cs@diagnostics$params = update_diagnostics(cs, step=i, leg="decay", out=output, runtime=a[1]+a[4], op=op)
    }
    #fit iota and rho given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic\n")
      a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_gauss_genomic_bam(
        biases = cs@biases, counts = cs@counts, design = cs@design,
        init = op$par, bf_per_kb = bf_per_kb, nthreads=ncores)))
      op=list(value=op.gen$value, par=c(op$par[c("beta_diag","beta_diag_centered","lambda_diag",
                                                 "log_decay","decay", "alpha")],
                                        op.gen$par[c("eC","eRJ","eDE",
                                                     "log_iota","log_rho","biases","lambda_iota","lambda_rho")]))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="bias", out=output, runtime=a[1]+a[4], op=op)
    }
    #make sure beta_diag is strictly increasing
    op$par$beta_diag = guarantee_beta_diag_increasing(op$par$beta_diag)
    if (fit.disp==T) {
      #fit exposures and dispersion
      if (verbose==T) cat("Gibbs",i,": Remaining parameters\n")
      a=system.time(output <- capture.output(op.disp <- csnorm:::csnorm_gauss_dispersion(
        biases = cs@biases, counts = subcounts, design = cs@design, init=op$par,
        dmin=dmin, dmax=dmax, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=iter,
        weight=subcounts.weight, init_alpha=init_alpha)))
      op=list(value=op.disp$value, par=c(op$par[c("beta_diag","beta_diag_centered","log_decay","decay",
                                                  "log_iota","log_rho","biases",
                                                  "lambda_iota","lambda_rho", "lambda_diag")],
                                         op.disp$par[c("eC","eRJ","eDE","alpha")]))
      cs@diagnostics$params = update_diagnostics(cs, step=i, leg="disp", out=output, runtime=a[1]+a[4], op=op)
      if (verbose==T) cat("Gibbs",i,": log-likelihood = ",op$value,"\n")
    }
  }
  if (verbose==T) cat("Done\n")
  op$par$init=init
  op$par$value=op$value
  cs@par=op$par
  cs@diagnostics$plot=ggplot(cs@diagnostics$params[,.(step,leg,value,out.last)])+
    geom_line(aes(step,value))+geom_point(aes(step,value,colour=out.last))+facet_wrap(~leg, scales = "free")
  return(cs)
}

