#' @include binless.R
NULL

#' Compute initial iota and rho biases assuming a poisson model
#' @keywords internal
#' 
initial_guess_decay = function(cs, cts.common, pseudocount=1e-2) {
  decay=cts.common[,.(log_decay=log((pseudocount+weighted.mean(count,weight))/(weighted.mean(exp(lmu.nosig),weight))),
                      ncounts=sum(weight)),keyby=c("name","dbin")]
  decay[,log_decay:=log_decay-weighted.mean(log_decay,ncounts),by=name]
  decay[,ncounts:=NULL]
  cs@par$decay=decay
  return(cs)
}

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
gauss_decay_muhat_mean = function(cs, cts.common) {
  cts.common[,kappaij:=eC+log_decay]
  dbins=cs@settings$dbins
  csd = cts.common[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
                      kappahat=weighted.mean(z+kappaij, weight/var),
                      std=1/sqrt(sum(weight/(2*var))), weight=sum(weight)/2), keyby=c("name", "dbin")] #each count appears twice
  return(csd)
}


#' Optimize decay parameters
#' @keywords internal
#' 
gauss_decay_optimize = function(csd, design, Kdiag, original_lambda_diag,
                                verbose=T, max_perf_iteration=1000, convergence_epsilon=1e-9) {
  Totalcbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  #cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  TotalDsets = design[,.N]
  XD=as.array(design[,decay])
  
  decay_out = list(beta_diag = matrix(0,nrow=TotalDsets,ncol=Kdiag),
                   beta_diag_diff = matrix(0,nrow=TotalDsets,ncol=Kdiag-2),
                   beta_diag_centered = c(), log_decay = c(), log_mean_counts = c(),
                   eC = c(),value = 0)
  decay_out$lambda_diag = c()
  for(uXD in unique(XD)) {
    if (verbose==T) cat("  group",uXD,"\n")
    Dsets = 0
    cbegin = c()
    for (d in 1:TotalDsets) {
      if(XD[d] == uXD) {
        cbegin = c(cbegin,Totalcbegin[d],Totalcbegin[d+1])
        Dsets = Dsets + 1
      }
    }
    SD = cbegin[2]-cbegin[1]
    cutsites = csd[,distance][cbegin[1]:(cbegin[2]-1)]
    X = as.matrix(generate_spline_base(log(cutsites), Kdiag))
    W=csd[,weight][cbegin[1]:(cbegin[2]-1)]
    diags = list(rep(1,Kdiag), rep(-2,Kdiag))
    D = bandSparse(Kdiag-2, Kdiag, k=0:2, diagonals=c(diags, diags[1]))
    if (!is.null(original_lambda_diag)) {
      lambda_diag = original_lambda_diag[uXD]
    } else {
      lambda_diag = 1
    }
    U_e=Matrix(rep.int(1,SD),ncol=1)
    kappa_hat=csd[,kappahat][cbegin[1]:(cbegin[2]-1)]
    sdl=csd[,std][cbegin[1]:(cbegin[2]-1)]
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        ncutsites = csd[,distance][cbegin[(2*d-1)]:(cbegin[2*d]-1)]
        nBsp = as.matrix(generate_spline_base(log(ncutsites), Kdiag))
        cutsites = c(cutsites,ncutsites)
        X = rbind(X,nBsp)
        SDd = cbegin[2*d]-cbegin[(2*d-1)]
        W = c(W,c(rep(0,SDd)))
        nU_e = Matrix(rep.int(1,SDd),ncol=1)
        U_e = bdiag(U_e,nU_e) 
        kappa_hat=c(kappa_hat,csd[,kappahat][cbegin[(2*d-1)]:(cbegin[2*d]-1)])
        sdl=c(sdl,csd[,std][cbegin[(2*d-1)]:(cbegin[2*d]-1)])
        SD = SD + SDd
      }
    } 
    
    S_m2 = Diagonal(x=1/sdl^2)
    
    Xt = cbind(U_e,X)
    Dt = cbind(matrix(0,ncol=Dsets, nrow=Kdiag-2),D)
    tmp_X_S_m2_X = crossprod(Diagonal(x=1/sdl)%*%Xt)
    tmp_X_S_m2_k = t(Xt)%*%S_m2%*%kappa_hat
    DtD = crossprod(Dt)
    diags = list(rep(1,Kdiag), rep(-2,Kdiag))
    C=-bandSparse(Kdiag, Kdiag-1, k=c(0,-1),diagonals=list(diags[[1]],-diags[[1]]))
    Ct=rbind(matrix(0,nrow=Dsets,ncol=Kdiag), cbind(crossprod(X,W),C))
    
    epsilon = Inf
    maxiter = 0
    betat = array(1,Kdiag+Dsets)
    
    while(epsilon > convergence_epsilon && maxiter < max_perf_iteration) {
      At = tmp_X_S_m2_X + Kdiag^2*lambda_diag^2*DtD
      #salt=Kdiag^2*lambda_diag^2*0.1*bdiag(Diagonal(Dsets)*0,Diagonal(Kdiag))
      #At.PD = At + salt
      At.PD = nearPD(At)$mat  
      fit = solve.QP(At.PD, tmp_X_S_m2_k, -Ct, meq = 1)
      nbetat = fit$solution
      epsilon = max(abs(betat-nbetat))
      betat=nbetat
      
      eC=as.array(betat[1:Dsets])
      beta=betat[(Dsets+1):(Kdiag+Dsets)]
      
      lambda_diag = (Kdiag - 2)/((Kdiag^2)*crossprod(D%*%beta)+1)
      lambda_diag = sqrt(as.numeric(lambda_diag))
      
      #cat("step ",maxiter," epsilon ",epsilon," lambda_diag ", lambda_diag,"\n")
      maxiter = maxiter+1
    }
    if (verbose==T) cat("   step",maxiter-1,": lambda_diag",lambda_diag,"\n")
    
    log_decay = X%*%beta
    log_mean_counts = log_decay
    SD = cbegin[2]-cbegin[1]
    log_mean_counts[1:SD] = log_mean_counts[1:SD] + rep(eC[1],SD)
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        SDd = cbegin[2*d]-cbegin[(2*d-1)]
        log_mean_counts[(SD+1):(SD+SDd)] = log_mean_counts[(SD+1):(SD+SDd)] + rep(eC[d],SDd)
        SD = SD + SDd
      }
    }
    for (d2 in 1:TotalDsets) {
      if(XD[d2] == uXD) {
        decay_out$beta_diag[d2,] = as.array(beta)
        decay_out$beta_diag_diff[d2,] = as.array(D%*%beta)
        decay_out$beta_diag_centered = c(decay_out$beta_diag_centered,as.array(beta))
      }
    }
    decay_out$log_mean_counts = c(decay_out$log_mean_counts,as.array(log_mean_counts))
    decay_out$log_decay = c(decay_out$log_decay,log_decay)
    decay_out$eC = c(decay_out$eC,eC)
    decay_out$lambda_diag = c(decay_out$lambda_diag,lambda_diag)
    
    decay_out$value = decay_out$value + sum(dnorm(kappa_hat, mean = as.array(decay_out$log_mean_counts), sd = as.array(sdl), log = TRUE))
    
  }
  decay_out$eC=as.array(decay_out$eC)
  decay_out$beta_diag=as.matrix(decay_out$beta_diag)
  decay_out$beta_diag_diff=as.matrix(decay_out$beta_diag_diff)
  #make decay data table, reused at next call
  dmat=csd[,.(name,dbin,distance,kappahat,std,ncounts=weight,kappa=decay_out$log_mean_counts,log_decay=decay_out$log_decay)]
  setkey(dmat,name,dbin)
  decay_out$decay=dmat
  decay_out[c("beta_diag","beta_diag_diff","beta_diag_centered","log_mean_counts")]=NULL
  return(decay_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
gauss_decay = function(cs, cts.common, verbose=T, update.eC=T) {
  if (verbose==T) cat(" Decay\n")
  csd = binless:::gauss_decay_muhat_mean(cs, cts.common)
  #run optimization
  op = binless:::gauss_decay_optimize(csd, cs@design, cs@settings$Kdiag, cs@par$lambda_diag,
                                      verbose=verbose, max_perf_iteration=cs@settings$iter,
                                      convergence_epsilon = cs@par$tol_decay)
  #restrict tolerance if needed
  precision = max(abs(op$log_decay - cs@par$log_decay))
  cs@par$tol_decay = min(cs@par$tol_decay, max(cs@settings$tol,precision/10))
  #update par slot
  if (update.eC==F) op$eC=NULL
  cs@par=modifyList(cs@par, op)
  if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
  return(cs)
}

