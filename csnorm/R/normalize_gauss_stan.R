#' @include csnorm.R
NULL

#' Single-cpu simplified initial guess
#' @keywords internal
#' 
csnorm_gauss_guess_stan = function(biases, counts, design, lambda, dmin, dmax, bf_per_kb=1, bf_per_decade=20,
                              iter=10000, dispersion=10, ...) {
  nBiases=design[,uniqueN(genomic)]
  init=list(log_iota=array(0,dim=biases[,.N]), log_rho=array(0,dim=biases[,.N]),
            log_decay=array(0,dim=counts[,.N]), eC=array(0,dim=design[,.N]),
            alpha=dispersion, lambda_iota=array(lambda,dim=nBiases), lambda_rho=array(lambda,dim=nBiases))
  op=csnorm_gauss_genomic_stan(biases, counts, design, init, bf_per_kb=bf_per_kb, iter=iter, estimate.lambdas=F, ...)
  #add diagonal decay inits
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  Decays=design[,uniqueN(decay)]
  beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
  op$par=c(list(beta_diag=beta_diag, lambda_diag=array(1,dim=Decays), log_decay=rep(0,counts[,.N]),
                alpha=dispersion, lambda_iota=init$lambda_iota, lambda_rho=init$lambda_rho),
           op$par[c("eC","eRJ","eDE","beta_iota","beta_rho","log_iota","log_rho")])
  return(op)
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
csnorm_gauss_genomic_stan = function(biases, counts, design, init, bf_per_kb=1, iter=10000,
                                verbose=T, estimate.lambdas=T, ...) {
  #add bias informations to counts
  csub=copy(counts)
  csub[,log_decay:=init$log_decay]
  bsub=biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
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
  cts=rbind(csub[,.(id=id1,R=contact.close, L=contact.far,  muR=exp(lmu.close), muL=exp(lmu.far),
                    etaL=eC + log_iota1, etaR=eC + log_rho1)],
            csub[,.(id=id1,R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_iota1, etaR=eC + log_rho1)],
            csub[,.(id=id2,R=contact.far,   L=contact.close,muR=exp(lmu.far),   muL=exp(lmu.close),
                    etaL=eC + log_iota2, etaR=eC + log_rho2)],
            csub[,.(id=id2,R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_iota2, etaR=eC + log_rho2)])
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
            alpha=init$alpha)
  op=optimize_stan_model(model=stanmodels$gauss_genomic, data=data, iter=iter, verbose=verbose, init=init, ...)
  #make nice output table
  bout=cbind(biases,as.data.table(op$par[c("log_iota","log_rho","log_mean_DL","log_mean_DR","log_mean_RJ")]))
  cout=cbind(cts,as.data.table(op$par[c("log_iota","log_rho","log_mean_cleft","log_mean_cright")]))
  cout=merge(cout,bout[,.(name,id,pos)],by="id")
  bout=rbind(bout[,.(cat="dangling L", name, id, pos, log_iota, log_rho, z=dangling.L/exp(log_mean_DL)-1, std=sqrt(1/exp(log_mean_DL)+1/init$alpha))],
             bout[,.(cat="dangling R", name, id, pos, log_iota, log_rho, z=dangling.R/exp(log_mean_DR)-1, std=sqrt(1/exp(log_mean_DL)+1/init$alpha))],
             bout[,.(cat="rejoined", name, id, pos, log_iota, log_rho, z=rejoined/exp(log_mean_RJ)-1, std=sqrt(1/exp(log_mean_RJ)+1/init$alpha))])
  cout=rbind(cout[,.(cat="contact L", name, id, pos, log_iota, log_rho, z=etaLhat-log_mean_cleft, std=sdL)],
             cout[,.(cat="contact R", name, id, pos, log_iota, log_rho, z=etaRhat-log_mean_cright, std=sdR)])
  bout=rbind(bout,cout)
  op$par$biases=bout
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

#' Run approximate gibbs sampler on with a single starting condition
#' @inheritParams run_gauss_stan
#' @param fit.decay,fit.genomic,fit.disp boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @keywords internal
#' @export
#' 
run_gauss_stan_gibbs = function(cs, init, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10,
                           ngibbs = 3, iter=100000, fit.decay=T, fit.genomic=T, fit.disp=T,
                           verbose=T, ncounts=100000, init_alpha=1e-5) {
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
  if (length(init)==1) {
    if (verbose==T) cat("Initial guess\n")
    init.a = system.time(init.output <- capture.output(init.op <- csnorm:::csnorm_gauss_guess_stan(
      biases = cs@biases, counts = cs@counts, design = cs@design, lambda=init[[1]], dmin=dmin, dmax=dmax,
      bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, iter = iter, dispersion=10, init_alpha=init_alpha)))
    #abort silently if initial guess went wrong
    if (length(grep("Line search failed",tail(init.output,1)))>0) {
      init.op$par$value=-.Machine$double.xmax
      cs@par=init.op$par
      cs@diagnostics$params = update_diagnostics(cs, step=0, leg="bias", out=init.output,
                                                 runtime=init.a[1]+init.a[4], op=init.op)
      return(cs)
    }
  } else {
    init.a = system.time(NULL)
    init.output = ""
    init.op = list(par=init,value=NA)
  }
  cs@diagnostics$params = update_diagnostics(cs, step=0, leg="bias", out=init.output,
                                             runtime=init.a[1]+init.a[4], op=init.op)
  op=init.op
  #make sure beta_diag is strictly increasing
  op$par$beta_diag = guarantee_beta_diag_increasing(op$par$beta_diag)
  #gibbs sampling
  for (i in 1:ngibbs) {
    #fit diagonal decay given iota and rho
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay\n")
      a=system.time(output <- capture.output(op.diag <- csnorm:::csnorm_gauss_decay(
        biases = cs@biases, counts = cs@counts, design=cs@design,
        init=op$par, dmin = dmin, dmax = dmax,
        bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter, init_alpha=init_alpha)))
      op=list(value=op.diag$value, par=c(op.diag$par[c("eC","beta_diag","beta_diag_centered",
                                                       "log_decay","decay", "lambda_diag")],
                                         op$par[c("eRJ","eDE","beta_iota","beta_rho", "alpha",
                                                  "lambda_iota","lambda_rho","log_iota","log_rho","biases")]))
      cs@diagnostics$params = update_diagnostics(cs, step=i, leg="decay", out=output, runtime=a[1]+a[4], op=op)
    }
    #fit iota and rho given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic\n")
      a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_gauss_genomic_stan(
        biases = cs@biases, counts = cs@counts, design = cs@design,
        init = op$par, bf_per_kb = bf_per_kb, iter = iter, init_alpha=init_alpha)))
      op=list(value=op.gen$value, par=c(op$par[c("beta_diag","beta_diag_centered","lambda_diag",
                                                 "log_decay","decay", "alpha")],
                                        op.gen$par[c("eC","eRJ","eDE","beta_iota","beta_rho",
                                                     "log_iota","log_rho","biases","lambda_iota","lambda_rho")]))
      cs@diagnostics$params = update_diagnostics(cs, step=i, leg="bias", out=output, runtime=a[1]+a[4], op=op)
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
                                                  "beta_iota","beta_rho","log_iota","log_rho","biases",
                                                  "lambda_iota","lambda_rho", "lambda_diag")],
                                         op.disp$par[c("eC","eRJ","eDE","alpha")]))
      cs@diagnostics$params = update_diagnostics(cs, step=i, leg="disp", out=output, runtime=a[1]+a[4], op=op)
      if (verbose==T) cat("Gibbs",i,": log-likelihood = ",op$value,"\n")
    }
  }
  if (verbose==T) cat("Done\n")
  op$par$init=init.op$par
  op$par$value=op$value
  cs@par=op$par
  cs@diagnostics$plot=ggplot(cs@diagnostics$params[,.(step,leg,value,out.last)])+
    geom_line(aes(step,value))+geom_point(aes(step,value,colour=out.last))+facet_wrap(~leg, scales = "free")
  return(cs)
}

#' Cut-site normalization (simplified gibbs sampling)
#' 
#' Alternates two approximations to the exact model, fitting the diagonal decay
#' and iota/rho.
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
#' @param lambdas positive numeric. Length scales to try out as initial
#'   condition.
#' @param ngibbs positive integer. Number of gibbs sampling iterations.
#' @param iter positive integer. Number of optimization steps for each stan 
#'   optimization call.
#' @param ncores positive integer. Number of cores to parallelize on.
#' @param verbose Display progress if TRUE
#' @param init_alpha positive numeric, default 1e-5. Initial step size of LBFGS
#'   line search.
#' @prefix character. If given, will save individual optimizations to files of
#'   form prefix_lambdaxxx.RData, where xxx is a float corresponding to the
#'   initial value of lambda. Useful in conjunction with \code{\link{recover_normalization}}
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
run_gauss_stan = function(cs, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10, lambdas=c(0.1,1,10),
                     ngibbs = 3, iter=100000, ncores=1, verbose=T, init_alpha=1e-5, prefix=NULL) {
  cs@binned=list() #erase old data if available
  cs@diagnostics=list()
  cs@par=list()
  registerDoParallel(cores=ncores)
  cs = foreach (lambda=lambdas, .combine=function(x,y){if (x@par$value[1]<y@par$value[1]){return(y)}else{return(x)}}) %dopar% {
    a=system.time(cs <- run_gauss_stan_gibbs(cs, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade,
                           bins_per_bf=bins_per_bf, init=lambda, ngibbs = ngibbs,
                           iter=iter, fit.decay=T, fit.genomic=T, fit.disp=T,
                           verbose=verbose, init_alpha=1e-5))
    cs@diagnostics$runtime=a[1]+a[4]
    if (!is.null(prefix)) save(cs, file=paste0(prefix,"_lambda",lambda,".RData"))
    cs
  }
  return(cs)
}

