#' @include binless.R
NULL

#' Get the first ncounts/d of each of d datasets, including zeros 
#' 
#' count is always a bit smaller because we censor those that are <dmin without adding more counts
#' 
#' @keywords internal
#' 
subsample_counts = function(cs, ncounts, dset=NA) {
  ncounts_per_dset=as.integer(ncounts/cs@design[,.N])
  if (is.na(dset)) {
    cts = foreach (d=cs@design[,name],.combine=rbind) %do% subsample_counts(cs,ncounts_per_dset,dset=d)
  } else {
    #get name and id of counts
    nbiases = cs@biases[name==dset, .N]
    ncounts = min(nbiases*(nbiases-1)/2,ncounts)
    ids=c(cs@biases[name==dset,.(minid=min(id),maxid=max(id))])
    cts=data.table(name=dset,id1=ids$minid,id2=(ids$minid+1):ids$maxid)
    while(cts[,.N]<ncounts)
      cts=rbind(cts,data.table(name=dset,id1=cts[.N,id1+1],id2=cts[.N,id1+2]:ids$maxid))
    cts=cts[1:ncounts]
    #merge positions and compute distances
    cts = merge(cts,cs@biases[,.(name,id,pos)],by.x=c("name","id1"),by.y=c("name","id"))
    cts = merge(cts,cs@biases[,.(name,id,pos)],by.x=c("name","id2"),by.y=c("name","id"),suffixes=c("1","2"))
    cts[,distance:=abs(pos2-pos1)]
    if (cs@settings$circularize>0) cts[,distance:=pmin(distance, cs@settings$circularize+1-distance)] 
    cts = cts[distance>=cs@settings$dmin]
    #merge counts
    setkey(cts, id1, id2, name)
    cts = merge(cts, cs@counts[,.(name,id1,id2,contact.close,contact.down,contact.far,contact.up)], all.x=T)
    cts[is.na(contact.close),contact.close:=0]
    cts[is.na(contact.down),contact.down:=0]
    cts[is.na(contact.far),contact.far:=0]
    cts[is.na(contact.up),contact.up:=0]
    if (cts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2)
      stop("dataset too sparse, please increase ncounts")
    return(cts)
  }
}

#' Compute means for a given counts matrix
#' @keywords internal
compute_means = function(cs, counts) {
  #compute background
  init=cs@par
  bsub=merge(init$biases[cat == "contact L",.(genomic.grp=group,pos,log_iota=eta)],
             init$biases[cat == "contact R",.(genomic.grp=group,pos,log_rho=eta)], by=c("genomic.grp","pos"))
  cpos=merge(cbind(cs@design[,.(name,decay.grp=decay,genomic.grp=genomic)],eC=init$eC), counts, by="name",all.x=F,all.y=T)
  setnames(bsub,c("genomic.grp","pos1","log_iota1","log_rho1"))
  cpos = merge(cpos, bsub, all.x=T, all.y=F, by=c("genomic.grp","pos1"))
  setnames(bsub,c("genomic.grp","pos2","log_iota2","log_rho2"))
  cpos = merge(cpos, bsub, all.x=T, all.y=F, by=c("genomic.grp","pos2"))
  cpos[,c("bin1","bin2","dbin"):=
         list(cut(pos1, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(pos2, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cpos=init$decay[,.(decay.grp=group,dbin,log_decay)][cpos,,on=c("decay.grp","dbin")]
  #compute signal
  if (cs@par$signal[,.N]>0 && length(cs@settings$sbins)>2) {
    signal = binless:::get_signal_matrix(cs, resolution = cs@settings$base.res, groups=cs@experiments[,.(name,groupname=name)])
    signal = rbind(signal[,.(name,bin1,bin2,phi)],signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
    cpos = signal[cpos,,on=c("name","bin1","bin2")]
  } else {
    cpos[,phi:=0]
  }
  #assemble
  cpos[,log_mu.base:=eC + log_decay + phi]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  
  cpos=cpos[,.(id1,id2,name,pos1,pos2,distance,contact.close,contact.far,contact.up,contact.down,
               log_decay,lmu.close,lmu.far,lmu.up,lmu.down)]
  setkeyv(cpos,key(cs@counts))
  cpos
}

#' Initial guess for normalization
#' @keywords internal
#' @export
#' 
optimize_stan_model = function(model, data, iter, verbose, init, ...) {
  out=capture.output(op<-optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...))
  cat(out,sep="\n")
  if (length(grep("Line search failed",tail(out,1)))>0) {
    op=optimizing(model, data=data, as_vector=F, hessian=F, iter=2, verbose=verbose, init=init, algorithm="Newton", ...)
    cat("!!! Line search error occurred, performed 2 steps of Newton optimization\n")
  }
  return(op)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
gauss_dispersion = function(cs, counts, weight=cs@design[,.(name,wt=1)], verbose=T) {
  if (verbose==T) cat(" Dispersion\n")
  #predict all means and put into table
  counts = binless:::compute_means(cs,counts)
  log_biases = dcast(cs@par$biases[cat%in%c("contact L","contact R"),.(group,pos,cat,eta)],group+pos~cat,value.var = "eta")
  setnames(log_biases,c("group","pos","log_iota","log_rho"))
  log_biases = merge(log_biases,cs@design[,.(name,group=genomic)],by=c("group"), allow.cartesian = T)
  log_biases = merge(log_biases,cs@biases,by=c("name","pos"),all.x=F,all.y=T)
  stopifnot(cs@biases[,.N]==log_biases[,.N])
  #
  #fit dispersion and exposures
  if (verbose==T) cat("  predict\n")
  bbegin=c(1,cs@biases[,.(name,row=.I)][name!=shift(name),row],cs@biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=cs@design[,.N], Biases=cs@design[,uniqueN(genomic)], Decays=cs@design[,uniqueN(decay)],
               XB=as.array(cs@design[,genomic]), XD=as.array(cs@design[,decay]),
               SD=cs@biases[,.N], bbegin=bbegin,
               cutsitesD=cs@biases[,pos], rejoined=cs@biases[,rejoined],
               danglingL=cs@biases[,dangling.L], danglingR=cs@biases[,dangling.R],
               N=counts[,.N], cbegin=cbegin,
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight[,wt]),
               log_iota=log_biases[,log_iota], log_rho=log_biases[,log_rho],
               log_mean_cclose=counts[,lmu.close], log_mean_cfar=counts[,lmu.far],
               log_mean_cup=counts[,lmu.up], log_mean_cdown=counts[,lmu.down])
  init=list(eC_sup=as.array(counts[,log(mean(contact.close/exp(lmu.close))),by=name][,V1]),
            eRJ=as.array(log_biases[,.(name,frac=rejoined/exp((log_iota+log_rho)/2))][,log(mean(frac)),by=name][,V1]),
            eDE=as.array(log_biases[,.(name,frac=(dangling.L/exp(log_iota)+dangling.R/exp(log_rho))/2)][
              ,log(mean(frac)),by=name][,V1]))
  init$mu=mean(exp(init$eC_sup[1]+counts[name==name[1],lmu.close]))
  init$alpha=max(0.001,1/(var(counts[name==name[1],contact.close]/init$mu)-1/init$mu))
  init$mu=NULL
  out=capture.output(op<-optimize_stan_model(model=binless:::stanmodels$gauss_dispersion, tol_param=cs@par$tol_disp,
                                             data=data, iter=cs@settings$iter, verbose=verbose, init=init,
                                             init_alpha=1e-9))
  #restrict tolerance if needed
  precision = max(abs(c(op$par["alpha"],recursive=T) - c(cs@par["alpha"],recursive=T)))
  cs@par$tol_disp = min(cs@par$tol_disp, max(cs@settings$tol, precision/10))
  #update parameters
  cs@par=modifyList(cs@par, op$par["alpha"])
  #cs@par$eC=cs@par$eC+op$par$eC_sup
  #
  #compute log-posterior
  Krow=cs@settings$Krow
  Kdiag=cs@settings$Kdiag
  cs@par$value = op$value + (Krow-2)/2*sum(log(cs@par$lambda_iota/exp(1))+log(cs@par$lambda_rho/exp(1))) +
    (Kdiag-2)/2*sum(log(cs@par$lambda_diag/exp(1)))
  if (verbose==T) {
    cat("  fit: dispersion",cs@par$alpha,"\n")
    cat("  log-likelihood = ",cs@par$value,"\n")
  }
  return(cs)
}
