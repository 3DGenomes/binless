#' @include csnorm.R
NULL

#' Compute decay etahat and weights using data as initial guess
#' @keywords internal
#' 
fill_parameters = function(cs, dispersion=1, lambda=1, fit.decay=T, fit.genomic=T, fit.disp=T) {
  init=list(alpha=dispersion)
  if (fit.decay==F) {
    Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    Decays=cs@design[,uniqueN(decay)]
    beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
    init=c(init,list(beta_diag=beta_diag, lambda_diag=array(lambda,dim=Decays), log_decay=rep(0,cs@counts[,.N])))
  }
  if (fit.genomic==F) {
    nBiases=cs@design[,uniqueN(genomic)]
    init=c(init,list(log_iota=array(0,dim=cs@biases[,.N]), log_rho=array(0,dim=cs@biases[,.N]),
                     lambda_iota=array(lambda,dim=nBiases), lambda_rho=array(lambda,dim=nBiases)))
  }
  if (fit.genomic==F || fit.decay==F || fit.disp==F) {
    init=c(init,list(eC=array(0,dim=cs@design[,.N]),  eRJ=array(0,dim=cs@design[,.N]),  eDE=array(0,dim=cs@design[,.N])))
  }
  cs@par=init
  cs
}

#' Compute decay etahat and weights using data as initial guess
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_data = function(cs, pseudocount=1e-2) {
  #bin distances
  stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
  dbins=10**seq(log10(cs@settings$dmin),log10(cs@settings$dmax)+stepsz,stepsz)
  csd=cs@counts[,.(name,distance,
                dbin=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12),
                kappahat=log(contact.close+contact.far+contact.up+contact.down+pseudocount),
                var=1/(contact.close+contact.far+contact.up+contact.down+pseudocount)+1/cs@par$alpha)]
  csd = csd[,.(distance=exp(mean(log(distance))), kappahat=sum(kappahat/var)/sum(1/var),
               std=1/sqrt(sum(1/var)), weight=4*.N), keyby=c("name", "dbin")]
  csd
}

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_mean = function(cs) {
  #add bias informations to counts
  init=cs@par
  csub=copy(cs@counts)
  csub[,log_decay:=init$log_decay]
  bsub=cs@biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  csub=merge(bsub[,.(id1=id,log_iota,log_rho)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_iota,log_rho)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(cs@design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
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
  stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
  dbins=10**seq(log10(cs@settings$dmin),log10(cs@settings$dmax)+stepsz,stepsz)
  csub[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #collect all counts in these bins
  csd = csub[,.(distance=exp(mean(log(distance))), kappahat=sum((z+kappaij)/var)/sum(1/var),
                std=1/sqrt(sum(1/var)), weight=4*.N), keyby=c("name", "dbin")]
  csd
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
csnorm_gauss_decay = function(cs, verbose=T, init.mean="mean", init_alpha=1e-7) {
  if (init.mean=="mean") {
    csd = csnorm:::csnorm_gauss_decay_muhat_mean(cs)
  } else {
    csd = csnorm_gauss_decay_muhat_data(cs)
  }
  #run optimization
  Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
  cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  data=list(Dsets=cs@design[,.N], Decays=cs@design[,uniqueN(decay)], XD=as.array(cs@design[,decay]),
            Kdiag=Kdiag, dmin=cs@settings$dmin, dmax=cs@settings$dmax, N=csd[,.N], cbegin=cbegin,
            kappa_hat=csd[,kappahat], sdl=csd[,std], dist=csd[,distance],
            weight=csd[,weight])
  #optimize from scratch, to avoid getting stuck. Slower but more robust
  op=optimize_stan_model(model=csnorm:::stanmodels$gauss_decay, data=data, iter=cs@settings$iter,
                         verbose=verbose, init=0, init_alpha=init_alpha)
  #make nice decay data table
  dmat=csd[,.(name,distance,kappahat,std,ncounts=weight,kappa=op$par$log_mean_counts)]
  setkey(dmat,name,distance)
  op$par$decay=dmat 
  #rewrite log_decay as if it were calculated for each count
  stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
  dbins=10**seq(log10(cs@settings$dmin),log10(cs@settings$dmax)+stepsz,stepsz)
  csub=cs@counts[,.(name,id1,id2,dbin=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  csd[,log_decay:=op$par$log_decay]
  a=csd[csub,.(name,id1,id2,log_decay),on=key(csd)]
  setkeyv(a,key(cs@counts))
  op$par$log_decay=a[,log_decay]
  #update par slot
  op$par$value=op$value
  op$par$beta_diag=guarantee_beta_diag_increasing(op$par$beta_diag)
  op$par$log_mean_counts=NULL
  cs@par=modifyList(cs@par, op$par)
  return(cs)
}

#' Compute genomic etahat and weights using data as initial guess
#' @keywords internal
#' 
csnorm_gauss_genomic_muhat_data = function(cs, pseudocount=1e-2) {
  dispersion=cs@par$alpha
  bts=rbind(cs@biases[,.(name,id,pos,cat="dangling L", etahat=log(dangling.L+pseudocount), std=sqrt(1/(dangling.L+pseudocount)+1/dispersion))],
            cs@biases[,.(name,id,pos,cat="dangling R", etahat=log(dangling.R+pseudocount), std=sqrt(1/(dangling.R+pseudocount)+1/dispersion))],
            cs@biases[,.(name,id,pos,cat="rejoined", etahat=log(rejoined+pseudocount), std=sqrt(1/(rejoined+pseudocount)+1/dispersion))])
  setkey(bts,id,name,cat)
  stopifnot(bts[,.N]==3*cs@biases[,.N])
  cts=rbind(cs@counts[,.(name,id=id1,pos=pos1, cat="contact R", etahat=log(contact.close+pseudocount), var=(1/(contact.close+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id1,pos=pos1, cat="contact L", etahat=log(contact.far+pseudocount), var=(1/(contact.far+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id1,pos=pos1, cat="contact R", etahat=log(contact.down+pseudocount), var=(1/(contact.down+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id1,pos=pos1, cat="contact L", etahat=log(contact.up+pseudocount), var=(1/(contact.up+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id2,pos=pos2, cat="contact R", etahat=log(contact.far+pseudocount), var=(1/(contact.far+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id2,pos=pos2, cat="contact L", etahat=log(contact.close+pseudocount), var=(1/(contact.close+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id2,pos=pos2, cat="contact R", etahat=log(contact.down+pseudocount), var=(1/(contact.down+pseudocount)+1/dispersion))],
            cs@counts[,.(name,id=id2,pos=pos2, cat="contact L", etahat=log(contact.up+pseudocount), var=(1/(contact.up+pseudocount)+1/dispersion))])
  cts=cts[,.(etahat=sum(etahat/var)/sum(1/var),std=sqrt(2/sum(1/var))),keyby=c("id","pos","name","cat")]
  setkeyv(cts,c("id","name","cat"))
  stopifnot(cts[,.N]==2*cs@biases[,.N])
  return(list(cts=cts,bts=bts))
}

#' Compute genomic etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_genomic_muhat_mean = function(cs) {
  #compute bias means
  init=cs@par
  bsub=copy(cs@biases)
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  bsub=merge(cbind(cs@design[,.(name)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
  bsub[,c("lmu.DL","lmu.DR","lmu.RJ"):=list(eDE+log_iota,eDE+log_rho,eRJ+(log_iota+log_rho)/2)]
  bts=rbind(bsub[,.(name,id,pos,cat="dangling L", etahat=dangling.L/exp(lmu.DL)-1+lmu.DL,
                    std=sqrt(1/exp(lmu.DL)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="dangling R", etahat=dangling.R/exp(lmu.DR)-1+lmu.DR,
                    std=sqrt(1/exp(lmu.DR)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="rejoined", etahat=rejoined/exp(lmu.RJ)-1+lmu.RJ,
                    std=sqrt(1/exp(lmu.RJ)+1/init$alpha))])
  setkey(bts,id,name,cat)
  stopifnot(bts[,.N]==3*cs@biases[,.N])
  bsub=bsub[,.(id,log_iota,log_rho)]
  #add bias informations to counts
  csub=copy(cs@counts)
  csub[,c("distance","log_decay"):=list(NULL,init$log_decay)]
  csub=merge(bsub[,.(id1=id,log_iota,log_rho)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_iota,log_rho)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(cs@design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  rm(bsub)
  #compute means
  csub[,lmu.base:=eC + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(lmu.base+log_iota1+log_rho2,
                                                           lmu.base+log_rho1 +log_rho2,
                                                           lmu.base+log_rho1 +log_iota2,
                                                           lmu.base+log_iota1+log_iota2)]
  csub[,lmu.base:=NULL]
  #collect all counts on left/right side
  cts=rbind(csub[,.(name, id=id1, pos=pos1, R=contact.close, L=contact.far,  muR=exp(lmu.close), muL=exp(lmu.far),
                    etaL=eC + log_iota1, etaR=eC + log_rho1)],
            csub[,.(name, id=id1, pos=pos1, R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_iota1, etaR=eC + log_rho1)],
            csub[,.(name, id=id2, pos=pos2, R=contact.far,   L=contact.close,muR=exp(lmu.far),   muL=exp(lmu.close),
                    etaL=eC + log_iota2, etaR=eC + log_rho2)],
            csub[,.(name, id=id2, pos=pos2, R=contact.down,  L=contact.up,   muR=exp(lmu.down),  muL=exp(lmu.up),
                    etaL=eC + log_iota2, etaR=eC + log_rho2)])
  rm(csub)
  cts[,c("varL","varR"):=list(1/muL+1/init$alpha,1/muR+1/init$alpha)]
  cts=rbind(cts[,.(cat="contact L", etahat=sum((L/muL-1+etaL)/varL)/sum(1/varL), std=sqrt(2/sum(1/varL))),by=c("name","id","pos")],
            cts[,.(cat="contact R", etahat=sum((R/muR-1+etaR)/varR)/sum(1/varR), std=sqrt(2/sum(1/varR))),by=c("name","id","pos")])
  setkey(cts,id,name,cat)
  stopifnot(cts[,.N]==2*cs@biases[,.N])
  return(list(bts=bts,cts=cts))
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
csnorm_gauss_genomic = function(cs, verbose=T, init.mean="mean", init_alpha=1e-7) {
  if (init.mean=="mean") {
    a = csnorm:::csnorm_gauss_genomic_muhat_mean(cs)
  } else {
    a = csnorm_gauss_genomic_muhat_data(cs)
  }
  bts=a$bts
  cts=a$cts
  #run optimization
  Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
  bbegin=c(1,cs@biases[,.(name,row=.I)][name!=shift(name),row],cs@biases[,.N+1])
  data=list(Dsets=cs@design[,.N], Biases=cs@design[,uniqueN(genomic)], XB=as.array(cs@design[,genomic]),
            Krow=Krow, SD=cs@biases[,.N], bbegin=bbegin, cutsitesD=cs@biases[,pos],
            eta_hat_RJ=bts[cat=="rejoined",etahat], sd_RJ=bts[cat=="rejoined",std],
            eta_hat_DL=bts[cat=="dangling L",etahat], sd_DL=bts[cat=="dangling L",std],
            eta_hat_DR=bts[cat=="dangling R",etahat], sd_DR=bts[cat=="dangling R",std],
            eta_hat_L=cts[cat=="contact L",etahat], sd_L=cts[cat=="contact L",std],
            eta_hat_R=cts[cat=="contact R",etahat], sd_R=cts[cat=="contact R",std])
  op=optimize_stan_model(model=csnorm:::stanmodels$gauss_genomic, data=data, iter=cs@settings$iter,
                         verbose=verbose, init=0, init_alpha=init_alpha)
  #make nice output table
  bout=rbind(bts[cat=="dangling L",.(cat, name, id, pos, etahat, std, eta=op$par$log_mean_DL)],
             bts[cat=="dangling R",.(cat, name, id, pos, etahat, std, eta=op$par$log_mean_DR)],
             bts[cat=="rejoined",.(cat, name, id, pos, etahat, std, eta=op$par$log_mean_RJ)])
  cout=rbind(cts[cat=="contact L",.(cat, name, id, pos, etahat, std, eta=op$par$log_mean_cleft)],
             cts[cat=="contact R",.(cat, name, id, pos, etahat, std, eta=op$par$log_mean_cright)])
  bout=rbind(bout,cout)
  setkey(bout, id, name, cat)
  op$par$biases=bout
  #update par slot
  op$par$value=op$value
  op$par$log_mean_DL=NULL
  op$par$log_mean_DR=NULL
  op$par$log_mean_RJ=NULL
  op$par$log_mean_cleft=NULL
  op$par$log_mean_cright=NULL
  cs@par=modifyList(cs@par, op$par)
  return(cs)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_gauss_dispersion = function(cs, counts, weight=cs@design[,.(name,wt=1)], verbose=T, init_alpha=1e-7) {
  Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))* cs@settings$bf_per_decade)
  bbegin=c(1,cs@biases[,.(name,row=.I)][name!=shift(name),row],cs@biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=cs@design[,.N], Biases=cs@design[,uniqueN(genomic)], Decays=cs@design[,uniqueN(decay)],
               XB=as.array(cs@design[,genomic]), XD=as.array(cs@design[,decay]),
               SD=cs@biases[,.N], bbegin=bbegin,
               cutsitesD=cs@biases[,pos], rejoined=cs@biases[,rejoined],
               danglingL=cs@biases[,dangling.L], danglingR=cs@biases[,dangling.R],
               Kdiag=Kdiag, dmin=cs@settings$dmin, dmax=cs@settings$dmax,
               N=counts[,.N], cbegin=cbegin,
               cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight[,wt]), log_iota=cs@par$log_iota, log_rho=cs@par$log_rho,
               beta_diag=cs@par$beta_diag)
  op=optimize_stan_model(model=csnorm:::stanmodels$gauss_dispersion, data=data, iter=cs@settings$iter,
                         verbose=verbose, init=cs@par, init_alpha=init_alpha)
  #update par slot
  op$par$value=op$value
  cs@par=modifyList(cs@par, op$par[c("eC","eRJ","eDE","alpha","value")])
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
#' @param fit.decay,fit.genomic,fit.disp boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @param verbose Display progress if TRUE
#' @param init_alpha positive numeric, default 1e-5. Initial step size of LBFGS
#'   line search.
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
#' 
run_gauss = function(cs, init=NULL, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10,
                     ngibbs = 3, iter=10000, fit.decay=T, fit.genomic=T, fit.disp=T,
                     verbose=T, ncounts=100000, init_alpha=1e-7) {
  #clean object if dirty
  cs@par=list() #in case we have a weird object
  cs@binned=list()
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings[c("circularize","dmin","dmax")],
                  list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
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
    if (is.data.table(cs@diagnostics$params)) laststep = cs@diagnostics$params[,max(step)] else laststep = 0
    init$beta_diag = guarantee_beta_diag_increasing(init$beta_diag)
    init.mean="mean"
    init$decay=NULL
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
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_genomic(cs, init.mean=init.mean, init_alpha=init_alpha)))
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

