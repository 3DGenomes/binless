#' @include binless.R
NULL

#' Compute initial iota and rho biases assuming a poisson model
#' @keywords internal
#' 
initial_guess_genomic = function(cs, cts.common, pseudocount=1e-2) {
  iotas=cts.common[cat=="contact L",log((pseudocount+weighted.mean(count,weight))/(weighted.mean(exp(lmu.nosig),weight))),keyby=c("name","id1")]
  iotas[,V1:=V1-mean(V1),by=name]
  cs@par$log_iota=as.array(iotas[,V1])
  rhos=cts.common[cat=="contact R",log((pseudocount+weighted.mean(count,weight))/(weighted.mean(exp(lmu.nosig),weight))),keyby=c("name","id1")]
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
  bts=rbind(bsub[,.(name,id,pos,cat="dangling L", etahat=dangling.L/exp(lmu.DL)-1+lmu.DL,
                    std=sqrt(1/exp(lmu.DL)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="dangling R", etahat=dangling.R/exp(lmu.DR)-1+lmu.DR,
                    std=sqrt(1/exp(lmu.DR)+1/init$alpha))],
            bsub[,.(name,id,pos,cat="rejoined", etahat=rejoined/exp(lmu.RJ)-1+lmu.RJ,
                    std=sqrt(1/exp(lmu.RJ)+1/init$alpha))])
  setkey(bts,id,name,cat)
  stopifnot(bts[,.N]==3*cs@biases[,.N])
  #counts
  cts.common[,etaij:=eC+log_bias]
  cts = cts.common[,.(etahat=weighted.mean(z+etaij, weight/var), std=1/sqrt(sum(weight/(2*var)))),
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
  
  genomic_out = list(beta_iota_diff = matrix(0, nrow = TotalDsets, ncol = Krow-2), beta_rho_diff = matrix(0, nrow = TotalDsets, ncol = Krow-2),
                     beta_iota = c(), beta_rho = c(), log_iota = c(), log_rho = c(),
                     log_mean_RJ = c(), log_mean_DL = c(), log_mean_DR = c(), log_mean_cleft  = c(), log_mean_cright = c(),
                     eC = c(), eDE = c(), eRJ = c(), lambda_iota = c(), lambda_rho = c(),
                     value = 0)
  
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
    Bsp = generate_cubic_spline(cutsites, Krow, sparse=T)
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
    U_e = cbind(c(rep.int(1,SD),rep.int(0,4*SD)),c(rep.int(0,SD),rep.int(1,2*SD),rep.int(0,2*SD)),c(rep.int(0,3*SD),rep.int(1,2*SD)))
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
        nBsp  = generate_cubic_spline(ncutsites, Krow, sparse=T)
        cutsites = c(cutsites,ncutsites)
        nX = rbind(cbind(nBsp/2,nBsp/2),bdiag(nBsp,nBsp),bdiag(nBsp,nBsp))
        SDd = bbegin[2*d]-bbegin[(2*d-1)]
        X = rbind(X,nX)
        Bsp = bdiag(Bsp,nBsp)
        W = rbind(W,Matrix(0,nrow=5*SDd,ncol=2*ncol(centering)))
        nU_e = cbind(c(rep.int(1,SDd),rep.int(0,4*SDd)),c(rep.int(0,SDd),rep.int(1,2*SDd),rep.int(0,2*SDd)),c(rep.int(0,3*SDd),rep.int(1,2*SDd)))
        U_e = bdiag(U_e,nU_e) 
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
      #
      tmp_Xt_Sm2_U = crossprod(X,S_m2%*%U_e)
      beta1 = solve(cholA, tmp_Xt_Sm2_U)
      beta2 = tmp_Lm1XtW %*% solve(cholWtXAm1XtW, crossprod(tmp_Xt_W,beta1))
      beta2 = solve(cholA,solve(cholA,beta2,system="DLt"),system="Pt") 
      beta_U = beta1-beta2
      #    
      e=solve(t(U_e)%*%S_m2%*%(U_e-X%*%beta_U),t(U_e)%*%S_m2%*%(etas-X%*%beta_y))
      
      nbeta = beta_y - beta_U%*%e
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
    
    eRJ = array(0,dim=c(Dsets))
    eDE = array(0,dim=c(Dsets))
    eC = array(0,dim=c(Dsets))
    for (d in 1:Dsets) {
      eRJ[d] = e[(3*(d-1)+1)]
      eDE[d] = e[(3*(d-1)+2)]
      eC[d] = e[(3*(d-1)+3)]
    }
    
    log_mean = U_e %*% e + X%*%beta
    SD = bbegin[2]-bbegin[1]
    log_iota = Bsp[1:SD,1:Krow]%*%beta_iota
    log_rho = Bsp[1:SD,1:Krow]%*%beta_rho
    
    log_mean_RJ = log_mean[1:SD]
    log_mean_DL = log_mean[(SD+1):(2*SD)]
    log_mean_DR = log_mean[(2*SD+1):(3*SD)]
    log_mean_cleft  = log_mean[(3*SD+1):(4*SD)]
    log_mean_cright = log_mean[(4*SD+1):(5*SD)]
    SD = 5*SD
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        SDd = bbegin[2*d]-bbegin[(2*d-1)]
        
        nlog_iota = Bsp[(SD/5+1):(SD/5+SDd),((d-1)*Krow+1):(d*Krow)]%*%beta_iota
        nlog_rho = Bsp[(SD/5+1):(SD/5+SDd),((d-1)*Krow+1):(d*Krow)]%*%beta_rho
        
        log_iota = c(as.array(log_iota),as.array(nlog_iota))
        log_rho = c(as.array(log_rho),as.array(nlog_rho))
        
        log_mean_RJ = c(log_mean_RJ,log_mean[(SD+1):(SD+SDd)])
        log_mean_DL = c(log_mean_DL,log_mean[(SD+(SDd+1)):(SD+(2*SDd))])
        log_mean_DR = c(log_mean_DR,log_mean[(SD+(2*SDd+1)):(SD+(3*SDd))])
        log_mean_cleft  = c(log_mean_cleft,log_mean[(SD+(3*SDd+1)):(SD+(4*SDd))])
        log_mean_cright = c(log_mean_cright,log_mean[(SD+(4*SDd+1)):(SD+(5*SDd))])
        
        SD = SD + 5*SDd
      }
    }
    
    genomic_out$beta_iota = c(genomic_out$beta_iota,as.array(beta_iota))
    genomic_out$beta_rho = c(genomic_out$beta_rho,as.array(beta_rho))
    for (d2 in 1:TotalDsets) {
      if(XB[d2] == uXB) {
        genomic_out$beta_iota_diff[d2,] = as.array(D1%*%beta_iota)
        genomic_out$beta_rho_diff[d2,] = as.array(D1%*%beta_rho)
      }
    }
    genomic_out$lambda_iota = c(genomic_out$lambda_iota,lambda_iota)
    genomic_out$lambda_rho = c(genomic_out$lambda_rho,lambda_rho)
    genomic_out$log_iota = c(genomic_out$log_iota,as.array(log_iota))
    genomic_out$log_rho = c(genomic_out$log_rho,as.array(log_rho))
    genomic_out$log_mean_RJ = c(genomic_out$log_mean_RJ,as.array(log_mean_RJ))
    genomic_out$log_mean_DL = c(genomic_out$log_mean_DL,as.array(log_mean_DL))
    genomic_out$log_mean_DR = c(genomic_out$log_mean_DR,as.array(log_mean_DR))
    genomic_out$log_mean_cleft  = c(genomic_out$log_mean_cleft,as.array(log_mean_cleft))
    genomic_out$log_mean_cright = c(genomic_out$log_mean_cright,as.array(log_mean_cright))
    genomic_out$eC = as.array(c(genomic_out$eC,eC))
    genomic_out$eRJ = as.array(c(genomic_out$eRJ,eRJ))
    genomic_out$eDE = as.array(c(genomic_out$eDE,eDE))
    genomic_out$beta_iota_diff = as.matrix(genomic_out$beta_iota_diff)
    genomic_out$beta_rho_diff = as.matrix(genomic_out$beta_rho_diff)
    
    genomic_out$value = genomic_out$value+sum(dnorm(etas, mean = as.array(log_mean), sd = as.array(sds), log = TRUE))
  }
  
  #make nice output table
  bout=rbind(bts[cat=="dangling L",.(cat, name, id, pos, etahat, std, eta=as.array(genomic_out$log_mean_DL))],
             bts[cat=="dangling R",.(cat, name, id, pos, etahat, std, eta=as.array(genomic_out$log_mean_DR))],
             bts[cat=="rejoined",.(cat, name, id, pos, etahat, std, eta=as.array(genomic_out$log_mean_RJ))])
  cout=rbind(cts[cat=="contact L",.(cat, name, id, pos, etahat, std, eta=as.array(genomic_out$log_mean_cleft))],
             cts[cat=="contact R",.(cat, name, id, pos, etahat, std, eta=as.array(genomic_out$log_mean_cright))])
  bout=rbind(bout,cout)
  setkey(bout, id, name, cat)
  genomic_out$biases=bout
  genomic_out[c("beta_iota_diff","beta_rho_diff","log_mean_RJ",
                "log_mean_DL","log_mean_DR","log_mean_cleft","log_mean_cright")]=NULL
  return(genomic_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
gauss_genomic = function(cs, cts.common, verbose=T, update.eC=T, constrain=F) {
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
  if (update.eC==F) {
    op$eC=NULL
  }
  op$eRJ=NULL
  op$eDE=NULL
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
  Bsp = binless:::generate_cubic_spline(cutsites, Krow, sparse=T)
  Dsets = cs@experiments[,.N]
  X = bdiag(lapply(1:Dsets,function(x){Bsp}))
  stopifnot(ncol(X) == length(cs@par$beta_iota), ncol(X) == length(cs@par$beta_rho))
  dt = data.table(name=rep(cs@experiments[,name],each=length(cutsites)), pos=rep(cutsites,Dsets),
                  log_iota = (X %*% cs@par$beta_iota)[,1], log_rho = (X %*% cs@par$beta_rho)[,1])
  return(dt)
}

