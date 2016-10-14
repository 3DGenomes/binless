#' @include csnorm.R
NULL

#' Single-cpu simplified initial guess
#' @keywords internal
#' 
csnorm_gauss_guess = function(biases, counts, design, lambda, dmin, dmax, bf_per_kb=1, bf_per_decade=5,
                              iter=10000, dispersion=10, ...) {
  nBiases=design[,uniqueN(genomic)]
  init=list(log_nu=array(0,dim=biases[,.N]), log_delta=array(0,dim=biases[,.N]),
            log_decay=array(0,dim=counts[,.N]), eC=array(0,dim=design[,.N]),
            alpha=dispersion, lambda_nu=array(lambda,dim=nBiases), lambda_delta=array(lambda,dim=nBiases))
  op=csnorm_gauss_genomic(biases, counts, design, init, bf_per_kb=bf_per_kb, iter=iter, ...)
  #add diagonal decay inits
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  op$par=c(list(beta_diag=beta_diag, lambda_diag=array(1,dim=Decays), log_decay=rep(0,counts[,.N]),
                alpha=dispersion, lambda_nu=init$lambda_nu, lambda_delta=init$lambda_delta),
           op$par[c("eC","eRJ","eDE","beta_nu","beta_delta","log_nu","log_delta")])
  return(op)
}

#' Single-cpu simplified fitting for nu and delta
#' @keywords internal
#' 
csnorm_gauss_decay = function(biases, counts, design, init, dmin, dmax,
                              bf_per_decade=5, bins_per_bf=10, iter=10000, verbose=T, ...) {
  for (n in biases[,unique(name)]) {
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #add bias informations to counts
  setkey(counts, id1, id2)
  csub=copy(counts)
  csub[,log_decay:=init$log_decay]
  bsub=biases[,.(id)]
  bsub[,c("log_nu","log_delta"):=list(init$log_nu,init$log_delta)]
  csub=merge(bsub[,.(id1=id,log_nu,log_delta)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_nu,log_delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  #add z-score and sd variables
  csub[,log_mu.base:=eC + log_nu1 + log_nu2 + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_delta1-log_delta2,
                                                           log_mu.base-log_delta1-log_delta2,
                                                           log_mu.base+log_delta2-log_delta1,
                                                           log_mu.base+log_delta1+log_delta2)]
  csub[,c("kappaij","log_mu.base"):=list(eC+log_decay,NULL)]
  csub=rbind(csub[,.(name,id1,id2,distance,kappaij,count=contact.far,mu=exp(lmu.far))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.down,mu=exp(lmu.down))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.up,mu=exp(lmu.up))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.close,mu=exp(lmu.close))])
  csub[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  #bin distances
  stepsz=1/(bins_per_bf*bf_per_decade)
  dbins=10**seq(log10(dmin),log10(dmax)+stepsz,stepsz)
  csub[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #collect all counts in these bins
  csd = csub[,.(mdist=exp(mean(log(distance))), kappahatl=sum((z+kappaij)/var)/sum(1/var),
                sdl=1/sqrt(sum(1/var)), weight=.N), keyby=c("name", "dbin")]
  #run optimization
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  data=list(Dsets=design[,.N], Decays=design[,uniqueN(decay)], XD=as.array(design[,decay]),
            Kdiag=Kdiag, dmin=dmin, dmax=dmax, N=csd[,.N], cbegin=cbegin,
            kappa_hat=csd[,kappahatl], sdl=csd[,sdl], dist=csd[,mdist],
            alpha=init$alpha, lambda_diag=init$lambda_diag, weight=csd[,weight])
  op=optimizing(stanmodels$gauss_decay, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init,
                init_alpha=1e-5, ...)
  #make nice decay data table
  csd[,log_decay:=op$par$log_decay]
  dmat=csd[,.(name, dist=mdist, log_decay, kappahatl, std=sdl)]
  dmat=merge(cbind(design[,.(name)],eC=init$eC), dmat, by="name",all.x=F,all.y=T)
  dmat[,z:=kappahatl-eC-log_decay]
  dmat=dmat[,.(name,dist,log_decay,z,std)]
  setkey(dmat,name,dist)
  op$par$decay=dmat 
  #rewrite log_decay as if it was calculated for each count
  csub=counts[,.(name,id1,id2,dbin=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5))]
  a=csd[csub,.(id1,id2,log_decay),on=key(csd)]
  setkeyv(a,key(counts))
  op$par$log_decay=a[,log_decay]
  return(op)
}

#' Single-cpu simplified fitting for nu and delta
#' @keywords internal
#' 
csnorm_gauss_genomic = function(biases, counts, design, init, bf_per_kb=1, iter=10000, verbose=T, ...) {
  for (n in biases[,unique(name)]) {
    stopifnot(counts[name==n,.N]==biases[name==n,.N*(.N-1)/2]) #needs to be zero-filled
  }
  #add bias informations to counts
  csub=copy(counts)
  csub[,log_decay:=init$log_decay]
  bsub=biases[,.(id)]
  bsub[,c("log_nu","log_delta"):=list(init$log_nu,init$log_delta)]
  csub=merge(bsub[,.(id1=id,log_nu,log_delta)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_nu,log_delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  #compute means
  csub[,lmu.base:=eC + log_nu1 + log_nu2 + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(lmu.base + log_delta1 - log_delta2,
                                                           lmu.base - log_delta1 - log_delta2,
                                                           lmu.base - log_delta1 + log_delta2,
                                                           lmu.base + log_delta1 + log_delta2)]
  csub[,lmu.base:=NULL]
  #collect all counts on left/right side
  cts=rbind(csub[,.(id=id1,R=contact.close, L=contact.far,  muR=exp(lmu.close), muL=exp(lmu.far),
                    etaL=eC + log_nu1 + log_delta1, etaR=eC + log_nu1 - log_delta1)],
            csub[,.(id=id1,R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_nu1 + log_delta1, etaR=eC + log_nu1 - log_delta1)],
            csub[,.(id=id2,R=contact.far,   L=contact.close,muR=exp(lmu.far),   muL=exp(lmu.close),
                    etaL=eC + log_nu2 + log_delta2, etaR=eC + log_nu2 - log_delta2)],
            csub[,.(id=id2,R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_nu2 + log_delta2, etaR=eC + log_nu2 - log_delta2)])
  cts[,c("varL","varR"):=list(1/muL+1/init$alpha,1/muR+1/init$alpha)]
  cts=cts[,.(etaLhat=sum((L/muL-1+etaL)/varL)/sum(1/varL),
             etaRhat=sum((R/muR-1+etaR)/varR)/sum(1/varR),
             sdL=sqrt(2)/sqrt(sum(1/varL)),
             sdR=sqrt(2)/sqrt(sum(1/varR))),keyby=id]
  stopifnot(cts[,.N]==biases[,.N])
  #run optimization
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)], XB=as.array(design[,genomic]),
            Krow=Krow, SD=biases[,.N], bbegin=bbegin,
            cutsitesD=biases[,pos], rejoined=biases[,rejoined],
            danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
            eta_hat_L=cts[,etaLhat], eta_hat_R=cts[,etaRhat], sd_L=cts[,sdL], sd_R=cts[,sdR],
            alpha=init$alpha, lambda_nu=init$lambda_nu, lambda_delta=init$lambda_delta)
  op=optimizing(stanmodels$gauss_genomic, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose,
                init=init, init_alpha=1e-5, ...)
  #make nice output table
  bout=cbind(biases,as.data.table(list(log_nu=op$par$log_nu,log_delta=op$par$log_delta)))
  bout=merge(cbind(design[,.(name)],eC=op$par$eC,eRJ=op$par$eRJ,eDE=op$par$eDE), bout, by="name",all.x=F,all.y=T)
  bout=merge(cts,bout,by="id",all=T)
  bout[,c("mu.DL","mu.DR","mu.RJ","lmu.L","lmu.R"):=list(
    exp(eDE+log_nu+log_delta),exp(eDE+log_nu-log_delta),exp(eRJ+log_nu),
    eC+log_nu+log_delta,eC+log_nu-log_delta)]
  bout=rbind(bout[,.(cat="dangling L", name, id, pos, log_nu, log_delta, z=dangling.L/mu.DL-1, std=sqrt(1/mu.DL+1/init$alpha))],
             bout[,.(cat="dangling R", name, id, pos, log_nu, log_delta, z=dangling.R/mu.DR-1, std=sqrt(1/mu.DR+1/init$alpha))],
             bout[,.(cat="rejoined", name, id, pos, log_nu, log_delta, z=rejoined/mu.RJ-1, std=sqrt(1/mu.RJ+1/init$alpha))],
             bout[,.(cat="contact L", name, id, pos, log_nu, log_delta, z=etaLhat-lmu.L, std=sdL)],
             bout[,.(cat="contact R", name, id, pos, log_nu, log_delta, z=etaRhat-lmu.R, std=sdR)])
  op$par$biases=bout
  op
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_gauss_dispersion = function(biases, counts, design, dmin, dmax, init,
                                   bf_per_kb = 1, bf_per_decade=5, iter = 10000,
                                   weight=array(1,dim=design[,.N]), verbose=T, ...) {
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
               weight=as.array(weight), beta_nu=init$beta_nu, beta_delta=init$beta_delta, beta_diag=init$beta_diag)
  op=optimizing(stanmodels$gauss_dispersion, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose,
                init=init, init_alpha=1e-5, ...)
  op$par$decay=data.table(dist=data$dist, decay=exp(op$par$log_decay), key="dist")
  return(op)
}

#' Run approximate gibbs sampler on with a single starting condition
#' @inheritParams run_simplified
#' @param fit.decay,fit.genomic boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @keywords internal
#' @export
#' 
run_gauss_gibbs = function(cs, init, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10,
                           ngibbs = 3, iter=100000, fit.decay=T, fit.genomic=T, fit.disp=T, verbose=T, ncounts=100000) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
                                    iter=iter, ngibbs=ngibbs))
  #fill counts matrix and take subset
  cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases, circularize=cs@settings$circularize)
  if (cs@counts[,.N] > ncounts*2+1) {
    subcounts = cs@counts[,.SD[(1:ncounts+as.integer(.N/2))],by=name]
  } else {
    subcounts=cs@counts
  }
  if (subcounts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2) stop("dataset too sparse, please increase ncounts")
  subcounts.weight=cs@counts[,.N]/subcounts[,.N]
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
  if (length(init)==1) {
    if (verbose==T) cat("Initial guess\n")
    init.a = system.time(init.output <- capture.output(init.op <- csnorm:::csnorm_gauss_guess(
      biases = cs@biases, counts = cs@counts, design = cs@design, lambda=init[[1]], dmin=dmin, dmax=dmax,
      bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, iter = iter, dispersion=10)))
    #abort silently if initial guess went wrong
    if (length(grep("Line search failed",tail(init.output,1)))>0) {
      init.op$par$value=-.Machine$double.xmax
      cs@par=init.op$par
      return(cs)
    }
  } else {
    init.a = system.time(NULL)
    init.output = ""
    init.op = list(par=init,value=NA)
  }
  cs@diagnostics=list(out.init=init.output, runtime.init=init.a[1]+init.a[4], op.init=init.op)
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
        init=op$par, dmin = dmin, dmax = dmax,
        bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter)))
      op=list(value=op.diag$value, par=c(op.diag$par[c("eC","beta_diag","beta_diag_centered",
                                                       "log_decay","decay")],
                                         op$par[c("eRJ","eDE","beta_nu","beta_delta", "alpha","lambda_diag",
                                                  "lambda_nu","lambda_delta","log_nu","log_delta","biases")]))
      cs@diagnostics[[paste0("out.decay",i)]]=output
      cs@diagnostics[[paste0("runtime.decay",i)]]=a[1]+a[4]
      cs@diagnostics[[paste0("op.decay",i)]]=op.diag
    }
    #fit nu and delta given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic\n")
      a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_gauss_genomic(
        biases = cs@biases, counts = cs@counts, design = cs@design,
        init = op$par, bf_per_kb = bf_per_kb, iter = iter)))
      op=list(value=op.gen$value, par=c(op$par[c("beta_diag","beta_diag_centered","lambda_diag",
                                                 "log_decay","decay", "alpha","lambda_nu","lambda_delta")],
                                        op.gen$par[c("eC","eRJ","eDE","beta_nu","beta_delta",
                                                     "log_nu","log_delta","biases")]))
      cs@diagnostics[[paste0("out.bias",i)]]=output
      cs@diagnostics[[paste0("runtime.bias",i)]]=a[1]+a[4]
      cs@diagnostics[[paste0("op.bias",i)]]=op.gen
    }
    #make sure beta_diag is strictly increasing
    for (d in 2:length(op$par$beta_diag)) {
      if (abs(op$par$beta_diag[d]-op$par$beta_diag[d-1])<10*.Machine$double.eps) {
        op$par$beta_diag[d:length(op$par$beta_diag)]=op$par$beta_diag[d:length(op$par$beta_diag)]+10*.Machine$double.eps
      }
    }
    if (fit.disp==T) {
      #fit exposures, lambdas and dispersion
      if (verbose==T) cat("Gibbs",i,": Remaining parameters\n")
      a=system.time(output <- capture.output(op.disp <- csnorm:::csnorm_gauss_dispersion(
        biases = cs@biases, counts = subcounts, design = cs@design, init=op$par,
        dmin=dmin, dmax=dmax, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=iter,
        weight=subcounts.weight)))
      op=list(value=op.disp$value, par=c(op$par[c("beta_diag","beta_diag_centered","log_decay","decay",
                                                  "beta_nu","beta_delta","log_nu","log_delta","biases")],
                                         op.disp$par[c("eC","eRJ","eDE","alpha",
                                                       "lambda_nu","lambda_delta", "lambda_diag")]))
      cs@diagnostics[[paste0("out.disp",i)]]=output
      cs@diagnostics[[paste0("runtime.disp",i)]]=a[1]+a[4]
      cs@diagnostics[[paste0("op.disp",i)]]=op.disp
    }
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
                    bins_per_bf=bins_per_bf, init=lambda, ngibbs = ngibbs,
                    iter=iter, fit.decay=T, fit.genomic=T, fit.disp=T, verbose=verbose)
  return(cs)
}

