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
  csd = cts.common[,.(eC=weighted.mean(z+eC, nobs/var), std=1/sqrt(sum(nobs/(2*var)))), keyby=c("name")] #each count appears twice
  cs@par$eC = as.array(csd[,eC])
  cs@par$value = sum(dnorm(0, mean = 0, sd = csd[,std], log = TRUE))
  if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
  return(cs)
}

