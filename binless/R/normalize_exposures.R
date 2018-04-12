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
  #process biases
  init=cs@par
  bsub=copy(cs@biases)
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  bsub=merge(cbind(cs@design[,.(name)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
  bsub[,c("lmu.DL","lmu.DR","lmu.RJ"):=list(eDE+log_iota,eDE+log_rho,eRJ+(log_iota+log_rho)/2)]
  #eDE
  bDE=rbind(bsub[,.(name,eDE,id,pos,cat="dangling L", z=dangling.L/exp(lmu.DL)-1,
                    var=1/exp(lmu.DL)+1/init$alpha)],
            bsub[,.(name,eDE,id,pos,cat="dangling R", z=dangling.R/exp(lmu.DR)-1,
                    var=1/exp(lmu.DR)+1/init$alpha)])
  dt.eDE = bDE[,.(eDE=weighted.mean(z+eDE, 1/var), std=1/sqrt(sum(1/var))), keyby=c("name")]
  #eRJ
  bRJ = bsub[,.(name,eRJ,id,pos,cat="rejoined", z=rejoined/exp(lmu.RJ)-1,
                    var=1/exp(lmu.RJ)+1/init$alpha)]
  dt.eRJ = bRJ[,.(eRJ=weighted.mean(z+eRJ, 1/var), std=1/sqrt(sum(1/var))), keyby=c("name")]
  #eC
  csd = cts.common[,.(eC=weighted.mean(z+eC, nobs/var), std=1/sqrt(sum(nobs/(2*var)))), keyby=c("name")] #each count appears twice
  #report
  cs@par$eDE = as.array(dt.eDE[,eDE])
  cs@par$eRJ = as.array(dt.eRJ[,eRJ])
  cs@par$eC = as.array(csd[,eC])
  cs@par$value = sum(dnorm(0, mean = 0, sd = c(dt.eDE[,std],dt.eRJ[,std],csd[,std]), log = TRUE))
  if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
  return(cs)
}

