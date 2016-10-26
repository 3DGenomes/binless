#' @include csnorm.R
NULL

#' Single-cpu simplified initial guess
#' @keywords internal
#' 
csnorm_gauss_guess_bam = function(biases, counts, design, dmin, dmax, bf_per_kb=1, bf_per_decade=20,
                              dispersion=10, nthreads=1, verbose=T) {
  init=list(log_iota=array(0,dim=biases[,.N]), log_rho=array(0,dim=biases[,.N]),
            log_decay=array(0,dim=counts[,.N]), eC=array(0,dim=design[,.N]),
            eRJ=array(0,dim=design[,.N]), eDE=array(0,dim=design[,.N]), alpha=dispersion)
  op = csnorm:::csnorm_gauss_genomic_bam(biases, counts, design, init, bf_per_kb=bf_per_kb, verbose=verbose, init.mean="data", nthreads=nthreads)
  #add diagonal decay inits
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  op$par=c(list(beta_diag=beta_diag, lambda_diag=array(1,dim=Decays), log_decay=rep(0,counts[,.N]),
                alpha=dispersion),
           op$par[c("eC","eRJ","eDE","log_iota","log_rho")])
  op
}
  

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
csnorm_gauss_decay = function(biases, counts, design, init, dmin, dmax,
                              bf_per_decade=20, bins_per_bf=10, iter=10000, verbose=T, ...) {
  #add bias informations to counts
  setkey(counts, id1, id2)
  csub=copy(counts)
  csub[,log_decay:=init$log_decay]
  bsub=biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  csub=merge(bsub[,.(id1=id,log_iota,log_rho)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_iota,log_rho)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  #add z-score and sd variables
  csub[,log_mu.base:=eC + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  csub[,c("kappaij","log_mu.base"):=list(eC+log_decay,NULL)]
  csub=rbind(csub[,.(name,id1,id2,distance,kappaij,count=contact.far,mu=exp(lmu.far))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.down,mu=exp(lmu.down))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.up,mu=exp(lmu.up))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.close,mu=exp(lmu.close))])
  csub[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  #bin distances
  stepsz=1/(bins_per_bf*bf_per_decade)
  dbins=10**seq(log10(dmin),log10(dmax)+stepsz,stepsz)
  csub[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #collect all counts in these bins
  csd = csub[,.(mdist=exp(mean(log(distance))), kappahatl=sum((z+kappaij)/var)/sum(1/var),
                sdl=1/sqrt(sum(1/var)), weight=.N), keyby=c("name", "dbin")]
  #run optimization
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  data=list(Dsets=design[,.N], Decays=design[,uniqueN(decay)], XD=as.array(design[,decay]),
            Kdiag=Kdiag, dmin=dmin, dmax=dmax, N=csd[,.N], cbegin=cbegin,
            kappa_hat=csd[,kappahatl], sdl=csd[,sdl], dist=csd[,mdist],
            alpha=init$alpha, weight=csd[,weight])
  #optimize from scratch, to avoid getting stuck. Slower but more robust
  op=optimize_stan_model(model=csnorm:::stanmodels$gauss_decay, data=data, iter=iter, verbose=verbose, init=0, ...)
  #make nice decay data table
  csd[,log_decay:=op$par$log_decay]
  dmat=csd[,.(name,dist=mdist,log_decay,z=kappahatl-op$par$log_mean_counts,std=sdl,ncounts=4*weight)]
  setkey(dmat,name,dist)
  op$par$decay=dmat 
  #rewrite log_decay as if it were calculated for each count
  csub=counts[,.(name,id1,id2,dbin=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  a=csd[csub,.(id1,id2,log_decay),on=key(csd)]
  setkeyv(a,key(counts))
  op$par$log_decay=a[,log_decay]
  return(op)
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
csnorm_gauss_genomic_bam = function(biases, counts, design, init, bf_per_kb=1, verbose=T, init.mean="mean", ...) {
  #compute bias means
  bsub=copy(biases)
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  bsub=merge(cbind(design[,.(name)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
  bsub[,c("lmu.DL","lmu.DR","lmu.RJ"):=list(eDE+log_iota,eDE+log_rho,eRJ+(log_iota+log_rho)/2)]
  bts=rbind(bsub[,.(name,id,pos,cat="dangling L", lmu=lmu.DL, etahat=dangling.L/exp(lmu.DL)-1+lmu.DL,
                    std=sqrt(1/exp(lmu.DL)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="dangling R", lmu=lmu.DR, etahat=dangling.R/exp(lmu.DR)-1+lmu.DR,
                    std=sqrt(1/exp(lmu.DR)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="rejoined", lmu=lmu.RJ, etahat=rejoined/exp(lmu.RJ)-1+lmu.RJ,
                    std=sqrt(1/exp(lmu.RJ)+1/init$alpha))])
  stopifnot(bts[,.N]==3*biases[,.N])
  bsub=bsub[,.(id,log_iota,log_rho)]
  #add bias informations to counts
  csub=copy(counts)
  csub[,log_decay:=init$log_decay]
  csub=merge(bsub[,.(id1=id,log_iota,log_rho)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_iota,log_rho)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  #compute means
  csub[,lmu.base:=eC + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(lmu.base+log_iota1+log_rho2,
                                                           lmu.base+log_rho1 +log_rho2,
                                                           lmu.base+log_rho1 +log_iota2,
                                                           lmu.base+log_iota1+log_iota2)]
  csub[,lmu.base:=NULL]
  #collect all counts on left/right side
  cts=rbind(csub[,.(name,id=id1, pos=pos1, R=contact.close, L=contact.far,  muR=exp(lmu.close), muL=exp(lmu.far),
                    etaL=eC + log_iota1, etaR=eC + log_rho1)],
            csub[,.(name,id=id1, pos=pos1, R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_iota1, etaR=eC + log_rho1)],
            csub[,.(name,id=id2, pos=pos2, R=contact.far,   L=contact.close,muR=exp(lmu.far),   muL=exp(lmu.close),
                    etaL=eC + log_iota2, etaR=eC + log_rho2)],
            csub[,.(name,id=id2, pos=pos2, R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_iota2, etaR=eC + log_rho2)])
  rm(csub)
  cts[,c("varL","varR"):=list(1/muL+1/init$alpha,1/muR+1/init$alpha)]
  cts=rbind(cts[,.(cat="contact L", lmu=sum(etaL/varL)/sum(1/varL), etahat=sum((L/muL-1+etaL)/varL)/sum(1/varL),
                   std=sqrt(2)/sqrt(sum(1/varL))),by=c("name","id","pos")],
            cts[,.(cat="contact R", lmu=sum(etaR/varR)/sum(1/varR), etahat=sum((R/muR-1+etaR)/varR)/sum(1/varR),
                   std=sqrt(2)/sqrt(sum(1/varR))),by=c("name","id","pos")])
  stopifnot(cts[,.N]==2*biases[,.N])
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
  cat("optimize with BAM") #necessary so logging goes smoothly
  data = foreach(gen=design[,uniqueN(genomic)], .combine=rbind) %do% {
    dsets=design[genomic==gen,name]
    sdata=data[name %in% dsets]
    if (init.mean=="mean") {
      sdata[,initval:=lmu]
    } else {
      sdata[,initval:=etahat]
    }
    #formula differs if there is more than one dataset
    if (length(dsets)>1) {
      sdata[,dset:=as.integer(factor(name))-1]
      fit=bam(formula = etahat ~ dset+is.dangling+is.rejoined-1
                         +s(pos, by=iota.coef,bs="ps", m=2, k=Krow) +s(pos, by=rho.coef,bs="ps", m=2, k=Krow),
              data=sdata, family=gaussian(), weight=sdata[,weight], scale=1,
              discrete=T, samfrac=0.1, mustart=sdata[,initval], ...)
    } else {
      fit=bam(formula = etahat ~ is.dangling+is.rejoined-1
              +s(pos, by=iota.coef,bs="ps", m=2, k=Krow) +s(pos, by=rho.coef,bs="ps", m=2, k=Krow),
              data=sdata, family=gaussian(), weight=sdata[,weight], scale=1,
              discrete=T, samfrac=0.1, mustart=sdata[,initval], ...)
    }
    sdata[,gam:=fit$fitted.values]
    sdata
  }
  #return values as if optimized by stan
  setkeyv(data,key(biases))
  eC=data[cat%in%c("contact L","contact R"), .(eC=mean(gam)),by=name]
  eDE=data[cat%in%c("dangling L","dangling R"), .(eDE=mean(gam)),by=name]
  eRJ=data[cat=="rejoined", .(eRJ=mean(gam)),by=name]
  log_iota=data[cat=="contact L",.(id,log_iota=gam-mean(gam)),by=name]
  log_rho=data[cat=="contact R",.(id,log_rho=gam-mean(gam)),by=name]
  data=data[merge(log_iota,log_rho,by=key(data))]
  op=list(value=-1, par=list(eC=eC[,eC], eDE=eDE[,eDE], eRJ=eRJ[,eRJ], log_iota=log_iota[,log_iota], log_rho=log_rho[,log_rho]))
  op$par$biases = data[,.(cat, name, id, pos, etahat, std, eta=lmu)]
  op
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_gauss_dispersion = function(biases, counts, design, dmin, dmax, init,
                                   bf_per_kb = 1, bf_per_decade = 20, iter = 10000,
                                   weight=design[,.(name,wt=1)], verbose=T, ...) {
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
               weight=as.array(weight[,wt]), log_iota=init$log_iota, log_rho=init$log_rho,
               beta_diag=init$beta_diag)
  op=optimize_stan_model(model=stanmodels$gauss_dispersion, data=data, iter=iter, verbose=verbose, init=init, ...)
  return(op)
}

#' Cut-site normalization (fast approximation, using mgcv BAM)
#' 
#' Alternates two approximations to the exact model, fitting the diagonal decay
#' and iota/rho.
#' 
#' @param cs CSnorm object as returned by \code{\link{merge_cs_norm_datasets}}
#' @param If provided, should correspond to a previous cs@par slot which will be used as a starting point
#' @param bf_per_kb positive numeric. Number of cubic spline basis functions per
#'   kilobase, for genomic bias estimates. Small values make the optimization 
#'   easy, but makes the genomic biases stiffer.
#' @param bf_per_decade positive numeric. Number of cubic spline basis functions
#'   per distance decade (in bases), for diagonal decay. Default parameter 
#'   should suffice.
#' @param bins_per_bf positive integer. Number of distance bins to split basis 
#'   functions into. Must be sufficiently small so that the diagonal decay is 
#'   approximately constant in that bin.
#' @param ngibbs positive integer. Number of gibbs sampling iterations.
#' @param iter positive integer. Number of optimization steps for each stan 
#'   optimization call.
#' @param fit.decay,fit.genomic,fit.disp boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @param verbose Display progress if TRUE
#' @param ncounts positive integer. Number of counts to use for dispersion estimation.
#' @param init_alpha positive numeric, default 1e-5. Initial step size of LBFGS
#'   line search (decay and dispersion steps).
#' @param ncores positive integer. Number of cores to parallelize the fit using bam. Does not scale well.
#'   
#' @return A csnorm object
#' @export
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
    init.a=system.time(init.output <- capture.output(op <- csnorm:::csnorm_gauss_guess_bam(
      biases = cs@biases, counts = cs@counts, design = cs@design, dmin=dmin, dmax=dmax,
      bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, nthreads=ncores)))
    init.output = "Flat init"
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
                                                  "log_iota","log_rho","biases")]))
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
                                        op.gen$par[c("eC","eRJ","eDE", "log_iota","log_rho","biases")]))
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
                                                  "log_iota","log_rho","biases","lambda_diag")],
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

