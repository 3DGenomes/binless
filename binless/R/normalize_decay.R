#' @include binless.R
NULL

#' Compute initial iota and rho biases assuming a poisson model
#' @keywords internal
#' 
initial_guess_decay = function(cs, cts.common, pseudocount=1e-2) {
  design=cs@design[,.(name,group=decay)]
  decay=cts.common[design][,.(log_decay=log(pseudocount+weighted.mean(count,nobs)/weighted.mean(exp(lmu.nosig),nobs)),
                      nobs=sum(nobs)),keyby=c("group","dbin")]
  decay[,log_decay:=log_decay-weighted.mean(log_decay,nobs),by=group]
  cs@par$decay=decay
  return(cs)
}

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
gauss_decay_muhat_mean = function(cs, cts.common) {
  dbins=cs@settings$dbins
  design=cs@design[,.(name,group=decay)]
  csd = cts.common[design][,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
                      kappahat=weighted.mean(z+log_decay, nobs/var),
                      std=1/sqrt(sum(nobs/(2*var))), nobs=sum(nobs)/2), keyby=c("group", "dbin")] #each count appears twice
  return(csd)
}


#' Optimize decay parameters
#' @keywords internal
#' 
gauss_decay_optimize = function(csd, Kdiag, original_lambda_diag,
                                verbose=T, max_perf_iteration=1000, convergence_epsilon=1e-9) {
  all_log_decay = c()
  decay_out = list(value = 0, lambda_diag = c())
  for(uXD in csd[,unique(group)]) {
    if (verbose==T) cat("  group",uXD,"\n")
    cutsites = csd[group==uXD,distance]
    lcs = log(cutsites)
    X = as.matrix(binless:::generate_spline_base(lcs,min(lcs),max(lcs), Kdiag))
    nobs_vec=csd[group==uXD,nobs]
    diags = list(rep(1,Kdiag), rep(-2,Kdiag))
    D = bandSparse(Kdiag-2, Kdiag, k=0:2, diagonals=c(diags, diags[1]))
    if (!is.null(original_lambda_diag)) {
      lambda_diag = original_lambda_diag[uXD]
    } else {
      lambda_diag = 1
    }
    kappa_hat=csd[group==uXD,kappahat]
    sdl=csd[group==uXD,std]
    
    S_m2 = Diagonal(x=1/sdl^2)
    
    Xt = X
    Dt = D
    tmp_X_S_m2_X = crossprod(Diagonal(x=1/sdl)%*%Xt)
    tmp_X_S_m2_k = t(Xt)%*%S_m2%*%kappa_hat
    DtD = crossprod(Dt)
    diags = list(rep(1,Kdiag), rep(-2,Kdiag))
    C=-bandSparse(Kdiag, Kdiag-1, k=c(0,-1),diagonals=list(diags[[1]],-diags[[1]]))
    Ct=C
    
    epsilon = Inf
    maxiter = 0
    betat = array(1,Kdiag)
    
    while(epsilon > convergence_epsilon && maxiter < max_perf_iteration) {
      At = tmp_X_S_m2_X + Kdiag^2*lambda_diag^2*DtD
      #salt=Kdiag^2*lambda_diag^2*0.1*bdiag(Diagonal(Dsets)*0,Diagonal(Kdiag))
      #At.PD = At + salt
      At.PD = nearPD(At)$mat  
      fit = solve.QP(At.PD, tmp_X_S_m2_k, -Ct)
      nbetat = fit$solution
      epsilon = max(abs(betat-nbetat))
      betat=nbetat
      
      beta=betat
      
      lambda_diag = (Kdiag - 2)/((Kdiag^2)*crossprod(D%*%beta)+1)
      lambda_diag = sqrt(as.numeric(lambda_diag))
      
      #cat("step ",maxiter," epsilon ",epsilon," lambda_diag ", lambda_diag,"\n")
      maxiter = maxiter+1
    }
    if (verbose==T) cat("   step",maxiter-1,": lambda_diag",lambda_diag,"\n")
    
    log_decay = as.array(X%*%beta)
    decay_out$value = decay_out$value + sum(dnorm(kappa_hat, mean = as.array(log_decay), sd = as.array(sdl), log = TRUE))
    
    avg.val = weighted.mean(log_decay,nobs_vec)
    log_decay = log_decay - rep(avg.val,length(log_decay))
    
    all_log_decay = c(all_log_decay,log_decay)
    decay_out$lambda_diag = c(decay_out$lambda_diag,lambda_diag)
    
  }
  #make decay data table, reused at next call
  csd[,log_decay:=all_log_decay]
  decay_out$decay=csd
  return(decay_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
gauss_decay = function(cs, cts.common, verbose=T) {
  if (verbose==T) cat(" Decay\n")
  csd = binless:::gauss_decay_muhat_mean(cs, cts.common)
  #run optimization
  op = binless:::gauss_decay_optimize(csd, cs@settings$Kdiag, cs@par$lambda_diag,
                                      verbose=verbose, max_perf_iteration=cs@settings$iter,
                                      convergence_epsilon = cs@par$tol_decay)
  #restrict tolerance if needed
  precision = max(abs(op$decay[,log_decay] - cs@par$decay[,log_decay]))
  cs@par$tol_decay = min(cs@par$tol_decay, max(cs@settings$tol,precision/10))
  #update par slot
  cs@par=modifyList(cs@par, op)
  if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
  return(cs)
}

