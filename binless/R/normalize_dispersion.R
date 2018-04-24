#' @include binless.R
NULL

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
gauss_dispersion = function(cs, cts.common, verbose=T) {
  if (verbose==T) cat(" Dispersion\n")
  #predict all means and put into table
  bts = binless:::gauss_common_muhat_mean_biases(cs)[,.(cat,count,mu,nobs)]
  cts = cts.common[,.(cat,count,mu,nobs=nobs/2)]
  data = rbind(bts,cts)
  alpha = MASS::theta.ml(data[,count], data[,mu], data[,sum(nobs)], data[,nobs], limit=10, eps=cs@par$tol_disp)
  #restrict tolerance if needed
  precision = max(abs(alpha - cs@par$alpha))
  cs@par$tol_disp = min(cs@par$tol_disp, max(cs@settings$tol, precision/10))
  #update parameters
  cs@par$alpha = alpha
  #
  #compute log-posterior
  Krow=cs@settings$Krow
  Kdiag=cs@settings$Kdiag
  cs@par$value = sum(data[,nobs]*(dnbinom(x=data[,count], size=alpha, mu=data[,mu],log=T)))
  if (verbose==T) {
    cat("  fit: alpha",cs@par$alpha,"\n")
    cat("  log-likelihood = ",cs@par$value,"\n")
  }
  return(cs)
}
