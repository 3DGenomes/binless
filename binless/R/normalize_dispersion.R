#' @include binless.R
NULL


#' This function is adapted from MASS::theta.ml (GPL-3)
#' It allows to specify the initial value of theta, and removes superfluous checks
#' @keywords internal
#' 
theta.ml = function (y, mu, weights, init.theta, limit = 25, eps = .Machine$double.eps^0.25, trace = FALSE) 
{
  score <- function(n, th, mu, y, w)
    sum(w * (digamma(th + y) - digamma(th) + log(th) + 1 - log(th + mu) - (y + th)/(mu + th)))
  info <- function(n, th, mu, y, w)
    sum(w * (-trigamma(th + y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu + th)^2))
  n = sum(weights)
  #t0 <- init.theta
  t0 <- n/sum(weights * (y/mu - 1)^2)
  it <- 0
  del <- 1
  if (trace) 
    message(sprintf("theta.ml: iter %d 'theta = %f'", it, 
                    signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, 
                                                   mu, y, weights))
    t0 <- t0 + del
    if (trace) 
      message("theta.ml: iter", it, " theta =", signif(t0))
  }
  t0
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
gauss_dispersion = function(cs, cts.common, verbose=T, alpha.min=0.01, ncores=1) {
  if (verbose==T) cat(" Dispersion\n")
  #predict all means and put into table
  bts = binless:::gauss_common_muhat_mean_biases(cs)[,.(cat,pos,count,mu,nobs)]
  cts = cts.common[,.(cat,pos=pos1,count,mu,nobs=nobs/2)]
  data = rbind(bts,cts)
  data[,bin:=cut(pos, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
  #take subset of rows
  bins=data[,unique(bin)]
  nrows=min(length(bins),cs@settings$nrows.dispersion)
  bins=bins[unique(as.integer(seq.int(1,length(bins),length.out=cs@settings$nrows.dispersion)))]
  data=data[bin%in%bins]
  #compute dispersion on each selected row
  registerDoParallel(cores=ncores)
  alphas = foreach (b=data[,unique(bin)],.combine=rbind) %dopar% {
    alpha=binless:::theta.ml(data[bin==b,count], data[bin==b,mu], data[bin==b,nobs], cs@par$alpha)
    data.table(bin=b,alpha=alpha)
  }
  stopImplicitCluster()
  #ggplot(alphas)+geom_point(aes(bin,alpha))+scale_y_log10()
  alpha = max(alphas[,median(alpha)],alpha.min)
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
