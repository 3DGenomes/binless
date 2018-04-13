#' @include binless.R
NULL

#' Compute initial iota and rho biases assuming a poisson model
#' @keywords internal
#' 
initial_guess_genomic = function(cs, cts.common, pseudocount=1e-2) {
  iotas=cts.common[cat=="contact L",log(pseudocount+weighted.mean(count,nobs)/weighted.mean(exp(lmu.nosig),nobs)),keyby=c("name","id1")]
  iotas[,V1:=V1-mean(V1),by=name]
  cs@par$log_iota=as.array(iotas[,V1])
  rhos=cts.common[cat=="contact R",log(pseudocount+weighted.mean(count,nobs)/weighted.mean(exp(lmu.nosig),nobs)),keyby=c("name","id1")]
  rhos[,V1:=V1-mean(V1),by=name]
  cs@par$log_rho=as.array(rhos[,V1])
  return(cs)
}



#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
gauss_genomic_muhat_mean = function(cs, cts.common) {
  #biases
  init=cs@par
  bsub=copy(cs@biases)
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  bsub=merge(cbind(cs@design[,.(name)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
  bsub[,c("lmu.DL","lmu.DR","lmu.RJ"):=list(eDE+log_iota,eDE+log_rho,eRJ+(log_iota+log_rho)/2)]
  bts=rbind(bsub[,.(name,id,pos,cat="dangling L", etahat=dangling.L/exp(lmu.DL)-1+log_iota,
                    std=sqrt(1/exp(lmu.DL)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="dangling R", etahat=dangling.R/exp(lmu.DR)-1+log_rho,
                    std=sqrt(1/exp(lmu.DR)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="rejoined", etahat=rejoined/exp(lmu.RJ)-1+(log_iota+log_rho)/2,
                    std=sqrt(1/exp(lmu.RJ)+1/init$alpha))])
  setkey(bts,id,name,cat)
  stopifnot(bts[,.N]==3*cs@biases[,.N])
  #counts
  cts = cts.common[,.(etahat=weighted.mean(z+log_bias, nobs/var), std=1/sqrt(sum(nobs/(2*var)))),
                   keyby=c("id1","name","cat")]
  cts = merge(cs@biases[,.(name,id,pos)],cts,by.x=c("name","id"),by.y=c("name","id1"))
  setkey(cts, id, name, cat)
  return(list(bts=bts,cts=cts))
}

gauss_genomic_optimize = function(bts, cts, biases, design, Krow, sbins,
                                  original_lambda_iota, original_lambda_rho, verbose=T,
                                  max_perf_iteration=1000, convergence_epsilon=1e-5, constrain=F) {
  XB = as.array(design[,genomic])
  
  #run optimization
  Totalbbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
  TotalDsets = design[,.N]
  
  genomic_out = list(log_iota = c(), log_rho = c(), lambda_iota = c(), lambda_rho = c(), value = 0)
  
  for(uXB in unique(XB)) {
    if (verbose==T) cat("  group",uXB,"\n")
    Dsets = 0
    bbegin = c()
    for (d in 1:TotalDsets) {
      if(XB[d] == uXB) {
        bbegin = c(bbegin,Totalbbegin[d],Totalbbegin[d+1])
        Dsets = Dsets + 1
      }
    }
    SD = bbegin[2]-bbegin[1]
    cutsites = biases[,pos][bbegin[1]:(bbegin[2]-1)]
    Bsp = generate_spline_base(cutsites, min(cutsites), max(cutsites), Krow)
    X = rbind(cbind(Bsp/2,Bsp/2),bdiag(Bsp,Bsp),bdiag(Bsp,Bsp))
    if (constrain == T) {
      sbinned=biases[bbegin[1]:(bbegin[2]-1), cut(pos, sbins, ordered_result=T,
                                                  right=F, include.lowest=T,dig.lab=12)]
      centering=Matrix(model.matrix(~ 0+sbinned))
      stopifnot(dim(centering)==c(SD,length(sbins)-1))
    } else {
      centering=Matrix(rep(1,SD),ncol=1)
    }
    W=cbind(rbind(Matrix(0,nrow=3*SD,ncol=ncol(centering)),centering,Matrix(0,nrow=SD,ncol=ncol(centering))),
            rbind(Matrix(0,nrow=4*SD,ncol=ncol(centering)),centering))
    diags = list(rep(1,Krow), rep(-2,Krow))
    D1 = bandSparse(Krow-2, Krow, k=0:2, diagonals=c(diags, diags[1]))
    if (!is.null(original_lambda_iota) && !is.null(original_lambda_rho)) {
      lambda_iota = original_lambda_iota[uXB]
      lambda_rho = original_lambda_rho[uXB]
    } else {
      lambda_iota = 1 
      lambda_rho = 1 
    }
    eta_hat_RJ=bts[cat=="rejoined",etahat][bbegin[1]:(bbegin[2]-1)]
    sd_RJ=bts[cat=="rejoined",std][bbegin[1]:(bbegin[2]-1)]
    eta_hat_DL=bts[cat=="dangling L",etahat][bbegin[1]:(bbegin[2]-1)]
    sd_DL=bts[cat=="dangling L",std][bbegin[1]:(bbegin[2]-1)]
    eta_hat_DR=bts[cat=="dangling R",etahat][bbegin[1]:(bbegin[2]-1)]
    sd_DR=bts[cat=="dangling R",std][bbegin[1]:(bbegin[2]-1)]
    eta_hat_L=cts[cat=="contact L",etahat][bbegin[1]:(bbegin[2]-1)]
    sd_L=cts[cat=="contact L",std][bbegin[1]:(bbegin[2]-1)]
    eta_hat_R=cts[cat=="contact R",etahat][bbegin[1]:(bbegin[2]-1)]
    sd_R=cts[cat=="contact R",std][bbegin[1]:(bbegin[2]-1)]
    
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        ncutsites = biases[,pos][bbegin[(2*d-1)]:(bbegin[2*d]-1)]
        nBsp  = generate_spline_base(ncutsites, min(ncutsites), max(ncutsites), Krow)
        cutsites = c(cutsites,ncutsites)
        nX = rbind(cbind(nBsp/2,nBsp/2),bdiag(nBsp,nBsp),bdiag(nBsp,nBsp))
        SDd = bbegin[2*d]-bbegin[(2*d-1)]
        X = rbind(X,nX)
        Bsp = bdiag(Bsp,nBsp)
        W = rbind(W,Matrix(0,nrow=5*SDd,ncol=2*ncol(centering)))
        eta_hat_RJ=c(eta_hat_RJ,bts[cat=="rejoined",etahat][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        sd_RJ=c(sd_RJ,bts[cat=="rejoined",std][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        eta_hat_DL=c(eta_hat_DL,bts[cat=="dangling L",etahat][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        sd_DL=c(sd_DL,bts[cat=="dangling L",std][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        eta_hat_DR=c(eta_hat_DR,bts[cat=="dangling R",etahat][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        sd_DR=c(sd_DR,bts[cat=="dangling R",std][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        eta_hat_L=c(eta_hat_L,cts[cat=="contact L",etahat][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        sd_L=c(sd_L,cts[cat=="contact L",std][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        eta_hat_R=c(eta_hat_R,cts[cat=="contact R",etahat][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        sd_R=c(sd_R,cts[cat=="contact R",std][bbegin[(2*d-1)]:(bbegin[2*d]-1)])
        SD = SD + SDd
      }
    } 
    
    SD = bbegin[2]-bbegin[1]
    etas = c(eta_hat_RJ[1:SD],eta_hat_DL[1:SD],eta_hat_DR[1:SD],eta_hat_L[1:SD],eta_hat_R[1:SD])
    sds = c((1/(sd_RJ[1:SD])),(1/(sd_DL[1:SD])),(1/(sd_DR[1:SD])),(1/(sd_L[1:SD])),(1/(sd_R[1:SD])))
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        SDd = bbegin[2*d]-bbegin[(2*d-1)]
        etas = c(etas,eta_hat_RJ[(SD+1):(SD+SDd)],eta_hat_DL[(SD+1):(SD+SDd)],eta_hat_DR[(SD+1):(SD+SDd)],eta_hat_L[(SD+1):(SD+SDd)],eta_hat_R[(SD+1):(SD+SDd)])
        sds = c(sds,(1/(sd_RJ[(SD+1):(SD+SDd)])),(1/(sd_DL[(SD+1):(SD+SDd)])),(1/(sd_DR[(SD+1):(SD+SDd)])),(1/(sd_L[(SD+1):(SD+SDd)])),(1/(sd_R[(SD+1):(SD+SDd)])))
        SD = SD + SDd
      }
    }
    
    S_m2 = Diagonal(x=sds^2)
    tmp_X_S_m2_X = crossprod(Diagonal(x=sds)%*%X)
    tmp_Xt_W = crossprod(X,W)
    
    epsilon = Inf
    maxiter = 0
    beta = Matrix(1,2*Krow)
    while(epsilon > convergence_epsilon && maxiter < max_perf_iteration) {
      D = bdiag(lambda_iota*D1,lambda_rho*D1)
      DtD = crossprod(D)
      if (maxiter==0) {
        cholA = Cholesky(tmp_X_S_m2_X + Krow^2*DtD + Diagonal(dim(DtD)[1])/(2*1**2), LDL=F, super=NA,Imult=1e-20) 
        stopifnot(!isLDL(cholA)) #do LLt cholesky so we can use crossprod for Gamma_v
      } else {
        cholA = update(cholA,tmp_X_S_m2_X + Krow^2*DtD + Diagonal(dim(DtD)[1])/(2*1**2))
      }
      tmp_Lm1XtW = solve(cholA,solve(cholA,tmp_Xt_W,system="P"),system="L") 
      if (maxiter==0) {
        cholWtXAm1XtW = Cholesky(crossprod(tmp_Lm1XtW), super=NA, Imult=1e-20)
      } else {
        cholWtXAm1XtW = update(cholWtXAm1XtW, t(tmp_Lm1XtW), mult=1e-20) #update providing M in A=MMt
      }
      #
      tmp_Xt_Sm2_etas = crossprod(X,S_m2%*%etas) #Kx1
      beta1 = solve(cholA, tmp_Xt_Sm2_etas)
      beta2 = tmp_Lm1XtW %*% solve(cholWtXAm1XtW, crossprod(tmp_Xt_W,beta1))
      beta2 = solve(cholA,solve(cholA,beta2,system="DLt"),system="Pt") 
      beta_y = beta1-beta2
      
      nbeta = beta_y
      epsilon = max(abs(beta-nbeta))
      beta = nbeta
      beta_iota = beta[1:Krow]
      beta_rho = beta[(Krow+1):(2*Krow)]
      
      lambda_iota = (Krow - 2)/((Krow**2)*crossprod(D1%*%beta_iota)+1e8)
      lambda_iota = sqrt(as.numeric(lambda_iota))
      lambda_rho = (Krow - 2)/((Krow**2)*crossprod(D1%*%beta_rho)+1e8)
      lambda_rho = sqrt(as.numeric(lambda_rho))
      
      maxiter = maxiter+1
    }
    if (verbose==T) cat("   step",maxiter-1,"lambda_iota",lambda_iota,"\n")
    
    log_mean = X%*%beta
    SD = bbegin[2]-bbegin[1]
    log_iota = Bsp[1:SD,1:Krow]%*%beta_iota
    log_rho = Bsp[1:SD,1:Krow]%*%beta_rho
    
    SD = 5*SD
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        SDd = bbegin[2*d]-bbegin[(2*d-1)]
        
        nlog_iota = Bsp[(SD/5+1):(SD/5+SDd),((d-1)*Krow+1):(d*Krow)]%*%beta_iota
        nlog_rho = Bsp[(SD/5+1):(SD/5+SDd),((d-1)*Krow+1):(d*Krow)]%*%beta_rho
        
        log_iota = c(as.array(log_iota),as.array(nlog_iota))
        log_rho = c(as.array(log_rho),as.array(nlog_rho))
        
        SD = SD + 5*SDd
      }
    }

    avg.iota = mean(log_iota)
    log_iota = log_iota - rep(avg.iota,length(log_iota))
    
    avg.rho = mean(log_rho)
    log_rho = log_rho - rep(avg.rho,length(log_rho))
    
    genomic_out$lambda_iota = c(genomic_out$lambda_iota,lambda_iota)
    genomic_out$lambda_rho = c(genomic_out$lambda_rho,lambda_rho)
    genomic_out$log_iota = c(genomic_out$log_iota,as.array(log_iota))
    genomic_out$log_rho = c(genomic_out$log_rho,as.array(log_rho))
    
    genomic_out$value = genomic_out$value+sum(dnorm(etas, mean = as.array(log_mean), sd = as.array(sds), log = TRUE))
  }
  
  #make nice output table
  bout=rbind(bts[cat=="dangling L",.(cat, name, id, pos, etahat, std, eta=genomic_out$log_iota)],
             bts[cat=="dangling R",.(cat, name, id, pos, etahat, std, eta=genomic_out$log_rho)],
             bts[cat=="rejoined",.(cat, name, id, pos, etahat, std, eta=(genomic_out$log_iota+genomic_out$log_rho)/2)])
  cout=rbind(cts[cat=="contact L",.(cat, name, id, pos, etahat, std, eta=genomic_out$log_iota)],
             cts[cat=="contact R",.(cat, name, id, pos, etahat, std, eta=genomic_out$log_rho)])
  bout=rbind(bout,cout)
  setkey(bout, id, name, cat)
  genomic_out$biases=bout
  return(genomic_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
gauss_genomic = function(cs, cts.common, verbose=T, constrain=F) {
  if (verbose==T) cat(" Genomic\n")
  a = binless:::gauss_genomic_muhat_mean(cs, cts.common)
  #run optimization
  op = binless:::gauss_genomic_optimize(a$bts, a$cts, cs@biases, cs@design, cs@settings$Krow, cs@settings$sbins,
                                        cs@par$lambda_iota, cs@par$lambda_rho, verbose=verbose,
                                        max_perf_iteration=cs@settings$iter,
                                        convergence_epsilon=cs@par$tol_genomic,
                                        constrain=constrain)
  #restrict tolerance if needed
  precision = max(max(abs(op$log_iota - cs@par$log_iota)), max(abs(op$log_rho - cs@par$log_rho)))
  cs@par$tol_genomic = min(cs@par$tol_genomic, max(cs@settings$tol, precision/10))
  #update par slot
  cs@par=modifyList(cs@par, op)
  if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
  return(cs)
}

#' Generate iota and rho genomic biases on evenly spaced points along the genome
#'
#' @param cs normalized CSnorm object
#' @param points_per_kb number of evaluation points per kb
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
generate_genomic_biases = function(cs, points_per_kb=10) {
  cutsites = cs@biases[, seq(min(pos), max(pos), by = 1000/points_per_kb)]
  Krow = cs@settings$Krow
  Bsp = binless:::generate_spline_base(cutsites, Krow)
  Dsets = cs@experiments[,.N]
  X = bdiag(lapply(1:Dsets,function(x){Bsp}))
  stopifnot(ncol(X) == length(cs@par$beta_iota), ncol(X) == length(cs@par$beta_rho))
  dt = data.table(name=rep(cs@experiments[,name],each=length(cutsites)), pos=rep(cutsites,Dsets),
                  log_iota = (X %*% cs@par$beta_iota)[,1], log_rho = (X %*% cs@par$beta_rho)[,1])
  return(dt)
}

