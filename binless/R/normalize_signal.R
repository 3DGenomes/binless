#' @include binless.R
NULL

#' Prepare for concurrent signal estimation 
#' @keywords internal
prepare_first_signal_estimation = function(biases, names, base.res) {
  ### build matrix
  #create an empty matrix containing all cells, even those with no cut-site intersection
  sbins=seq(biases[,min(pos)-1],biases[,max(pos)+1+base.res],base.res)
  signal.bins=unique(cut(c(sbins,head(sbins,n=-1)+base.res/2), sbins,
                         ordered_result=T, right=F, include.lowest=T,dig.lab=12))
  signal.mat=CJ(name=names,bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
  signal.mat[,c("phi","beta"):=list(0,0)]
  setkey(signal.mat,name,bin1,bin2)
  stopifnot(all(signal.mat[,.N,by=name]$N==signal.mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  return(list(signal=signal.mat,sbins=sbins))
}

#' fit signal using sparse fused lasso
#' @keywords internal
#' 
gauss_signal_muhat_mean = function(cs, cts.common) {
  cts = cts.common[,.(name,bin1=pmin(bin1,bin2),bin2=pmax(bin1,bin2),count,z,var,mu,lmu.nosig,
                      log_decay,log_bias,nobs=nobs/2)]
  stopifnot(cts[,all(bin1<=bin2)])
  return(cts)
}

#get metadata for signal calculation
#consists of outliers, which should be discarded for signal detection,
#and information relative to constraining the signal wrt the decay
#' @keywords internal
#' 
get_signal_metadata = function(cs, cts, resolution) {
  cts_compressed = cts[,.(nobs=sum(nobs), z=weighted.mean(z,nobs/var),var=1/weighted.mean(1/var,nobs)),
                       keyby=c("name","bin1","bin2")]
  #biases
  bad.biases=rbind(cts_compressed[,.(bin1,bin2,nobs,var,z)],cts_compressed[,.(bin1=bin2,bin2=bin1,nobs,var,z)])[
    ,.(z=sum(nobs*z/var)/sum(nobs^2/var^2)),by=bin1]
  bad.biases[,z:=scale(z)]
  bad.biases[,is.out:=-abs(z)<qnorm(cs@settings$qmin)]
  #ggplot(bad.biases)+geom_point(aes(bin1,z,colour=is.out))+geom_hline(aes(yintercept=qnorm(cs@settings$qmin)))
  bad.rows=bad.biases[is.out==T,bin1]
  if (bad.biases[,sum(is.out)/.N>0.1]) cat(" Warning: removing ",bad.biases[,100*sum(is.out)/.N],"% of all rows!\n")
  #decay
  cts_compressed[,diag.idx:=unclass(bin2)-unclass(bin1)]
  # bad.decays=cts_compressed[,.(z=sum(nobs*z/var)/sum(nobs^2/var^2)),by=diag.idx]
  # bad.decays[,z:=scale(z)]
  # bad.decays[,is.out:=-abs(z)<qnorm(cs@settings$qmin)]
  diag.rm = ceiling(cs@settings$dmin/resolution)
  # bad.decays[diag.idx<=diag.rm,is.out:=T]
  # #ggplot(bad.decays)+geom_point(aes(diag.idx,z,colour=is.out))
  # bad.diagonals=bad.decays[is.out==T,diag.idx]
  # if (bad.decays[,sum(is.out)/.N>0.1]) cat(" Warning: removing ",bad.decays[,100*sum(is.out)/.N],"% of all counter diagonals!\n")
  #orthogonality
  dx = 1.01*(log10(cs@settings$dmax)-log10(cs@settings$dmin))/(cs@settings$Kdiag-3)
  orth=data.table(diag.idx=0:cts_compressed[,max(diag.idx)],
                  segment=ceiling((log10(0:cts_compressed[,max(diag.idx)]*resolution)-log10(cs@settings$dmin))/(dx/2))) #double constraint
  orth[,rank:=frank(segment,ties.method = "dense")-1]
  return(list(bad.diagonals=0:diag.rm, bad.rows=bad.rows, diag.grp=orth[,rank]))
  #return(list(bad.diagonals=bad.diagonals,bad.rows=bad.rows, diag.grp=orth[,rank]))
}


#' make a summarized cts matrix used for one IRLS iteration
#' @keywords internal
#' 
compress_cts = function(cs, cts) {
  cts.new = cts[,.(phihat=weighted.mean(count/mu-1,nobs/(2*var)), nobs=sum(nobs/(2*var)),
                   log_decay=weighted.mean(log_decay,nobs)), keyby=c("name","bin1","bin2")]
  cts.new[,count:=phihat+2] #or any strictly positive number (but beware of exp(lmu.nosig) )
  cts.new = cs@par$signal[,.(name,bin1,bin2,phi)][cts.new]
  cts.new[,.(name,bin1,bin2,count,lmu.nosig=log(count/(phihat+1))-phi, log_decay, nobs=2*nobs*(1/cs@par$alpha + (1+phihat)/count))]
}

#' fit signal using sparse fused lasso
#' @keywords internal
#' 
gauss_signal = function(cs, cts.common, verbose=T, ncores=1, fix.lambda1=F, fix.lambda1.at=NA,
                        fix.lambda2=F, fix.lambda2.at=NA) {
  if (verbose==T) cat(" Signal\n")
  #cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
  cts = binless:::gauss_signal_muhat_mean(cs, cts.common)
  metadata = binless:::get_signal_metadata(cs, cts, cs@settings$base.res)
  cts = binless:::compress_cts(cs, cts)
  #
  if (verbose==T) cat("  predict\n")
  #perform fused lasso on signal
  groupnames=cts[,unique(name)]
  nbins=length(cs@settings$sbins)-1
  csigs = foreach(g=groupnames) %do% {
    csig=new("CSbsig", mat=cs@par$signal[name==g], cts=cts[name==g],
             settings=list(metadata=metadata,
                           nbins=nbins, dispersion=cs@par$alpha,
                           last.beta=cs@par$signal[name==g,beta],
                           tol.val=cs@par$tol_signal, nperf=1)) #only one IRLS iteration
    csig@state = binless:::gfl_compute_initial_state(csig, diff=F)
    csig
  }
  registerDoParallel(cores=min(ncores,length(groupnames)))
  params = foreach(csig=csigs, .combine=rbind, .export=c("verbose","fix.lambda1","fix.lambda1.at",
                                                         "fix.lambda2","fix.lambda2.at")) %dopar% {
                                                           binless:::fused_lasso(csig, positive=T, fixed=F, constrained=F, verbose=verbose,
                                                                                 fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                                                                                 fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at)
                                                         }
  stopImplicitCluster()
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  #store new signal in cs and update eC
  #ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+
  #  scale_fill_gradient2()+coord_fixed()
  #ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=phi==0))
  setkey(mat,name,bin1,bin2)
  #restrict tolerance if needed
  precision = max(abs(mat[,beta]-cs@par$beta.phi))
  cs@par$tol_signal = min(cs@par$tol_signal, max(cs@settings$tol, precision/10))
  #set new parameters
  cs@par$signal=mat[,.(name,bin1,bin2,phihat,weight,ncounts,phi,beta,diag.grp,diag.idx)]
  cs@par$beta.phi=mat[,beta]
  params=merge(cbind(cs@design[,.(name)],eC=cs@par$eC), params, by="name",all=T)
  cs@par$eC=as.array(params[,eC+eCprime])
  cs@par$eCprime=as.array(params[,eCprime])
  cs@par$lambda1=as.array(params[,lambda1])
  cs@par$lambda2=as.array(params[,lambda2])
  cs@par$value = params[,sum(BIC)]
  if (verbose==T) {
    cat("  fit: lambda1",cs@par$lambda1[1],"lambda2",cs@par$lambda2[1],"\n")
    cat("  BIC = ",cs@par$value, "\n")
  }
  return(cs)
}
