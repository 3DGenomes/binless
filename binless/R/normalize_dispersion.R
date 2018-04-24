#' @include binless.R
NULL

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
