#' @include binless.R
NULL

#' Compute initial exposures assuming a poisson model
#' @keywords internal
#' 
initial_guess_exposures = function(cs, cts.common, pseudocount=1e-2) {
  #biases
  cs@par$eDE=as.array(cs@biases[,log(pseudocount+mean(dangling.L+dangling.R)/2),keyby=c("name")]$V1)
  cs@par$eRJ=as.array(cs@biases[,log(pseudocount+mean(rejoined)),keyby=c("name")]$V1)
  cs@par$eC=array(0,cs@experiments[,.N])
  cs@par$eC = as.array(cts.common[,log(pseudocount+weighted.mean(count,nobs)),keyby=c("name")]$V1)
  return(cs)
}

#' Compute eC using previous mean
#' @keywords internal
#' 
gauss_exposures = function(cs, cts.common, verbose=T) {
  if (verbose==T) cat(" Exposures\n")
  #biases
  init=cs@par
  bsub=copy(cs@biases)
  bsub=merge(cbind(cs@design[,.(name,group=genomic)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
  bts=rbind(bsub[,.(group,name,cat="rejoined",pos, count=rejoined,expo=eRJ,nobs=1)],
            bsub[,.(group,name,cat="dangling L",pos, count=dangling.L,expo=eDE,nobs=1)],
            bsub[,.(group,name,cat="dangling R",pos, count=dangling.R,expo=eDE,nobs=1)])
  setkey(bts,group,cat,pos)
  bts = bts[init$biases[cat%in%c("rejoined","dangling L","dangling R"),.(group,cat,pos,eta)]]
  bts[,mu:=exp(expo+eta)]
  bts[,c("z","var"):=list(count/mu-1,var=1/mu+1/init$alpha)]
  #eDE
  dt.eDE = bts[cat!="rejoined",.(eDE=weighted.mean(z+expo, 1/var), std=1/sqrt(sum(1/var))), keyby=c("name")]
  #eRJ
  dt.eRJ = bts[cat=="rejoined",.(eRJ=weighted.mean(z+expo, 1/var), std=1/sqrt(sum(1/var))), keyby=c("name")]
  #eC
  dt.eC = cts.common[,.(eC=weighted.mean(z+eC, nobs/var), std=1/sqrt(sum(nobs/(2*var)))), keyby=c("name")] #each count appears twice
  #report
  cs@par$eDE = as.array(dt.eDE[,eDE])
  cs@par$eRJ = as.array(dt.eRJ[,eRJ])
  cs@par$eC = as.array(dt.eC[,eC])
  cs@par$value = sum(dnorm(0, mean = 0, sd = c(dt.eDE[,std],dt.eRJ[,std],dt.eC[,std]), log = TRUE))
  if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
  return(cs)
}

