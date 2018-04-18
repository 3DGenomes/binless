#' @include binless.R
NULL

#' Compute initial iota and rho biases assuming a poisson model
#' @keywords internal
#' 
initial_guess_genomic = function(cs, cts.common, pseudocount=1e-2) {
  design=cs@design[,.(name,group=genomic)]
  iotas=cts.common[design][cat=="contact L",.(log_iota=log(pseudocount+weighted.mean(count,nobs)/weighted.mean(exp(lmu.nosig),nobs))),keyby=c("group","pos1")]
  iotas[,log_iota:=log_iota-mean(log_iota),by=group]
  rhos=cts.common[design][cat=="contact R",.(log_rho=log(pseudocount+weighted.mean(count,nobs)/weighted.mean(exp(lmu.nosig),nobs))),keyby=c("group","pos1")]
  rhos[,log_rho:=log_rho-mean(log_rho),by=group]
  both = merge(iotas,rhos,all=T,by=c("group","pos1"))
  #
  biases = rbind(
    both[,.(cat="dangling L",group,pos=pos1,eta=log_iota)],
    both[,.(cat="dangling R",group,pos=pos1,eta=log_rho)],
    both[,.(cat="rejoined",group,pos=pos1,eta=(log_iota+log_rho)/2)],
    both[,.(cat="contact L",group,pos=pos1,eta=log_iota)],
    both[,.(cat="contact R",group,pos=pos1,eta=log_iota)]
  )
  setkey(biases,group,cat,pos)
  cs@par$biases=biases
  return(cs)
}



#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
gauss_genomic_muhat_mean = function(cs, cts.common) {
  #biases
  init=cs@par
  bsub=copy(cs@biases)
  bsub=merge(cbind(cs@design[,.(name,group=genomic)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
  bts=rbind(bsub[,.(group,cat="rejoined",pos, count=rejoined,expo=eRJ,nobs=1)],
            bsub[,.(group,cat="dangling L",pos, count=dangling.L,expo=eDE,nobs=1)],
            bsub[,.(group,cat="dangling R",pos, count=dangling.R,expo=eDE,nobs=1)])
  setkey(bts,group,cat,pos)
  bts = bts[init$biases[cat%in%c("rejoined","dangling L","dangling R"),.(group,cat,pos,eta)]]
  bts[,mu:=exp(expo+eta)]
  bts[,c("etahat","var"):=list(count/mu-1+eta,var=1/mu+1/init$alpha)]
  bts[,c("eta","mu","count","expo"):=NULL]
  #counts
  cts = cs@design[,.(name,group=genomic)][cts.common][,.(group,cat,pos=pos1,etahat=z+log_bias,var,nobs)]
  biasmat = rbind(bts,cts)
  biasmat[,cat:=ordered(cat,levels=c("rejoined","dangling L", "dangling R", "contact L", "contact R"))] #important to match design matrix X
  biasmat = biasmat[,.(etahat=weighted.mean(etahat, nobs/var), std=1/sqrt(sum(nobs/(2*var))), nobs=sum(nobs)),
                   keyby=c("group","cat","pos")]
  stopifnot(biasmat[,!any(is.na(cat))])
  stopifnot(all(biasmat[,.SD[,.N],by=pos]$V1==5))
  return(biasmat)
}

gauss_genomic_optimize = function(biasmat, design, Krow, sbins,
                                  original_lambda_iota, original_lambda_rho, verbose=T,
                                  max_perf_iteration=1000, convergence_epsilon=1e-5, constrain=F) {
  XB = as.array(design[,genomic])
  
  genomic_out = list(lambda_iota = c(), lambda_rho = c(), value = 0, genomic_beta = data.table(), biases = data.table())
  
  for(uXB in unique(XB)) {
    if (verbose==T) cat("  group",uXB,"\n")
    cutsites = biasmat[group==uXB&cat==cat[1],pos]
    Bsp = generate_spline_base(cutsites, min(cutsites), max(cutsites), Krow)
    X = rbind(cbind(Bsp/2,Bsp/2),bdiag(Bsp,Bsp),bdiag(Bsp,Bsp)) #because the biasmat is properly sorted and all 5 cats per pos are there
    SD=length(cutsites)
    if (constrain == T) {
      sbinned=cut(cutsites, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)
      centering=Matrix(model.matrix(~ 0+sbinned))
      stopifnot(dim(centering)==c(SD,length(sbins)-1))
    } else {
      centering=Matrix(biasmat[group==uXB,.(nobs=sum(nobs)),keyby=pos]$nobs,ncol=1)
      centering=centering/mean(centering)
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
    etas = biasmat[group==uXB,etahat]
    sds = biasmat[group==uXB,std]
    
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
    
    beta_iota = beta[1:Krow]
    beta_rho = beta[(Krow+1):(2*Krow)]
    genomic_out$genomic_beta = rbind(genomic_out$genomic_beta,
                                     data.table(group=uXB, bf=1:Krow, beta_iota=beta_iota, beta_rho=beta_rho))
    
    log_iota = as.array(Bsp%*%beta_iota)
    log_rho = as.array(Bsp%*%beta_rho)
    
    avg.iota = mean(log_iota)
    log_iota = log_iota - rep(avg.iota,length(log_iota))
    
    avg.rho = mean(log_rho)
    log_rho = log_rho - rep(avg.rho,length(log_rho))
    
    genomic_out$biases = rbind(genomic_out$biases,
                               biasmat[group==uXB,.(group,cat,pos,etahat,std,eta=c((log_iota+log_rho)/2,log_iota,log_rho,log_iota,log_rho),nobs)])
    
    genomic_out$lambda_iota = c(genomic_out$lambda_iota,lambda_iota)
    genomic_out$lambda_rho = c(genomic_out$lambda_rho,lambda_rho)
    
    genomic_out$value = genomic_out$value+sum(dnorm(etas, mean = as.array(X %*% beta), sd = as.array(sds), log = TRUE))
  }
  
  return(genomic_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
gauss_genomic = function(cs, cts.common, verbose=T, constrain=F) {
  if (verbose==T) cat(" Genomic\n")
  biasmat = binless:::gauss_genomic_muhat_mean(cs, cts.common)
  #run optimization
  op = binless:::gauss_genomic_optimize(biasmat, cs@design, cs@settings$Krow, cs@settings$sbins,
                                        cs@par$lambda_iota, cs@par$lambda_rho, verbose=verbose,
                                        max_perf_iteration=cs@settings$iter,
                                        convergence_epsilon=cs@par$tol_genomic,
                                        constrain=constrain)
  #restrict tolerance if needed
  precision = max(abs(op$biases[,eta] - cs@par$biases[,eta]))
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
  beta_iota = cs@par$genomic_beta[,beta_iota]
  beta_rho = cs@par$genomic_beta[,beta_rho]
  cutsites = cs@biases[, seq(min(pos), max(pos), by = 1000/points_per_kb)]
  Krow = cs@settings$Krow
  Bsp = binless:::generate_spline_base(cutsites, Krow)
  Ngroups = cs@design[,uniqueN(genomic)]
  X = bdiag(lapply(1:cs@design[,uniqueN(genomic)],function(x){Bsp}))
  stopifnot(ncol(X) == length(beta_iota), ncol(X) == length(beta_rho))
  dt = data.table(group=rep(1:Ngroups,each=length(cutsites)), pos=rep(cutsites,Ngroups),
                  log_iota = (X %*% beta_iota)[,1], log_rho = (X %*% beta_rho)[,1])
  return(dt)
}

