#' @include csnorm.R
NULL

#' Compute decay etahat and weights using data as initial guess
#' @keywords internal
#' 
fill_parameters_perf = function(cs, dispersion=10, lambda=1, fit.decay=T, fit.genomic=T, fit.disp=T) {
  nBiases=cs@design[,uniqueN(genomic)]
  Decays=cs@design[,uniqueN(decay)]
  init=list(alpha=dispersion)
  if (fit.decay==F) {
    Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
    init=c(init,list(beta_diag=beta_diag, lambda_diag=array(lambda,dim=Decays), log_decay=rep(0,cs@counts[,.N])))
  }
  if (fit.genomic==F) {
    init=c(init,list(log_iota=array(0,dim=cs@biases[,.N]), log_rho=array(0,dim=cs@biases[,.N]),
                     lambda_iota=array(lambda,dim=nBiases), lambda_rho=array(lambda,dim=nBiases)))
  }
  if (fit.genomic==F || fit.decay==F || fit.disp==F) {
    init=c(init,list(eC=array(0,dim=cs@design[,.N]),  eRJ=array(0,dim=cs@design[,.N]),  eDE=array(0,dim=cs@design[,.N])))
  }
  cs@par=init
  cs
}

#' Compute decay etahat and weights using data as initial guess
#' @keywords internal
#' 
fill_parameters_outer = function(cs, dispersion=10, lambda=1, fit.decay=T, fit.genomic=T, fit.disp=T) {
  nBiases=cs@design[,uniqueN(genomic)]
  Decays=cs@design[,uniqueN(decay)]
  init=list(alpha=dispersion,
            lambda_iota=array(lambda,dim=nBiases), lambda_rho=array(lambda,dim=nBiases), lambda_diag=array(lambda,dim=Decays))
  if (fit.decay==F) {
    Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
    init=c(init,list(beta_diag=beta_diag, log_decay=rep(0,cs@counts[,.N])))
  }
  if (fit.genomic==F) {
    init=c(init,list(log_iota=array(0,dim=cs@biases[,.N]), log_rho=array(0,dim=cs@biases[,.N])))
  }
  if (fit.genomic==F || fit.decay==F || fit.disp==F) {
    init=c(init,list(eC=array(0,dim=cs@design[,.N]),  eRJ=array(0,dim=cs@design[,.N]),  eDE=array(0,dim=cs@design[,.N])))
  }
  cs@par=init
  cs
}

#' Compute decay etahat and weights using data as initial guess
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_data = function(cs, zdecay, pseudocount=1e-2) {
  #bin distances
  dbins=cs@settings$dbins
  mcounts=melt(cs@counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0,.(name,distance,category,count)]
  mcounts[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #add zero counts
  mcounts = rbind(mcounts[,.(name,dbin,category,count,weight=1)], zdecay[,.(name,dbin=bdist,category=NA,count=0,weight=nzero)])
  #compute z-scores and sum counts
  mcounts[,c("kappahat","var"):=list(log(count+pseudocount), 1/(count+pseudocount)+1/cs@par$alpha)]
  csd = mcounts[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
                   kappahat=weighted.mean(kappahat, weight/var),
                   std=1/sqrt(sum(weight/var)), weight=sum(weight)), keyby=c("name", "dbin")]
  stopifnot(csd[,!is.na(distance)])
  return(csd)
}

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_mean = function(cs, zdecay) {
  #add bias informations to counts
  init=cs@par
  csub=copy(cs@counts)
  stopifnot(length(init$log_decay)==csub[,.N])
  csub[,log_decay:=init$log_decay]
  bsub=cs@biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  csub=merge(bsub[,.(id1=id,log_iota,log_rho)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_iota,log_rho)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(cs@design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  #add z-score and sd variables
  csub[,log_mu.base:=eC + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  csub[,c("kappaij","log_mu.base"):=list(eC+log_decay,NULL)]
  csub[,c("kappaij"):=list(eC+log_decay)]
  csub=rbind(csub[,.(name,id1,id2,distance,kappaij,count=contact.far,mu=exp(lmu.far))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.down,mu=exp(lmu.down))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.up,mu=exp(lmu.up))],
             csub[,.(name,id1,id2,distance,kappaij,count=contact.close,mu=exp(lmu.close))])
  csub[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  #bin distances
  dbins=cs@settings$dbins
  csub[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #collect all counts in these bins and add zeros
  stopifnot(!is.null(init$decay))
  csd = rbind(csub[count>0,.(name,dbin,z,kappaij,var,weight=1)],
              zdecay[init$decay,.(name,dbin=bdist,z=-1,kappaij=kappa,var=1/exp(kappa)+1/init$alpha,weight=nzero)])
  csd = csd[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
               kappahat=weighted.mean(z+kappaij, weight/var),
               std=1/sqrt(sum(weight/var)), weight=sum(weight)), keyby=c("name", "dbin")]
  stopifnot(csd[,!is.na(distance)])
  return(csd)
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
csnorm_gauss_decay = function(cs, zdecay, verbose=T, init.mean="mean", init_alpha=1e-7, type=c("outer","perf"), fit_model=c('stan','nostan'), max_perf_iteration=1000, convergence_epsilon=1e-5) {
  type=match.arg(type)
  fit_model=match.arg(fit_model)
  if (init.mean=="mean") {
    csd = csnorm:::csnorm_gauss_decay_muhat_mean(cs, zdecay)
  } else {
    csd = csnorm:::csnorm_gauss_decay_muhat_data(cs, zdecay)
  }
  #run optimization
  Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
  Totalcbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  #cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  TotalDsets = cs@design[,.N]
  XD=as.array(cs@design[,decay])
  
  all_beta_diag = matrix(0,nrow=TotalDsets,ncol=Kdiag)
  all_beta_diag_diff = matrix(0,nrow=TotalDsets,ncol=Kdiag-2)
  all_beta_diag_centered = c()
  all_log_decay = c()
  all_log_mean_counts = c()
  all_eC = c()
  all_value = 0
  if (type=="perf") {
    all_lambda_diag = c()
  }
  for(uXD in unique(XD)) {
    Dsets = 0
    cbegin = c()
    cbegin_stan = c(1)
    stanXD = c()
    for (d in 1:TotalDsets) {
      if(XD[d] == uXD) {
        cbegin = c(cbegin,Totalcbegin[d],Totalcbegin[d+1])
        cbegin_stan = c(cbegin_stan,(Totalcbegin[d+1]-Totalcbegin[d]+tail(cbegin_stan, n=1)))
        Dsets = Dsets + 1
        stanXD = c(stanXD,1) 
      }
    }
    SD = cbegin[2]-cbegin[1]
    #sparse spline design
    cutsites = csd[,distance][cbegin[1]:(cbegin[2]-1)]
    q = 3
    dx = 1.01*(max(log(cutsites))-min(log(cutsites)))/(Kdiag-q)
    t = min(log(cutsites)) - dx*0.01 + dx * seq(-q,Kdiag-q+3)
    X = spline.des(log(cutsites), knots = t, outer.ok = T, sparse=F)$design
    W=csd[,weight][cbegin[1]:(cbegin[2]-1)]
    stan_W = W
    diags = list(rep(1,Kdiag), rep(-2,Kdiag))
    D = bandSparse(Kdiag-2, Kdiag, k=0:2, diagonals=c(diags, diags[1]))
    if (type=="outer" || (!is.null(cs@par$lambda_diag))) {
      lambda_diag = cs@par$lambda_diag[uXD]
    } else {
      lambda_diag = 1
    }
    U_e=Matrix(rep.int(1,SD),ncol=1)
    kappa_hat=csd[,kappahat][cbegin[1]:(cbegin[2]-1)]
    sdl=csd[,std][cbegin[1]:(cbegin[2]-1)]
    if (Dsets > 1) {
      for (d in 2:Dsets) {
        ncutsites = csd[,distance][cbegin[(2*d-1)]:(cbegin[2*d]-1)]
        dx = 1.01*(max(log(ncutsites))-min(log(ncutsites)))/(Kdiag-q)
        t = min(log(ncutsites)) - dx*0.01 + dx * seq(-q,Kdiag-q+3)
        nBsp = spline.des(log(ncutsites), knots = t, outer.ok = T, sparse=F)$design
        cutsites = c(cutsites,ncutsites)
        X = rbind(X,nBsp)
        SDd = cbegin[2*d]-cbegin[(2*d-1)]
        W = c(W,c(rep(0,SDd)))
        stan_W=c(stan_W,csd[,weight][cbegin[(2*d-1)]:(cbegin[2*d]-1)])
        nU_e = Matrix(rep.int(1,SDd),ncol=1)
        U_e = bdiag(U_e,nU_e) 
        kappa_hat=c(kappa_hat,csd[,kappahat][cbegin[(2*d-1)]:(cbegin[2*d]-1)])
        sdl=c(sdl,csd[,std][cbegin[(2*d-1)]:(cbegin[2*d]-1)])
        SD = SD + SDd
      }
    } 
    if(fit_model=='stan') {
      
      data=list(Dsets=Dsets, Decays=1, XD=as.array(stanXD),
                Kdiag=Kdiag, dmin=cs@settings$dmin, dmax=cs@settings$dmax, N=SD, cbegin=as.array(cbegin_stan),
                kappa_hat=kappa_hat, sdl=sdl, dist=cutsites,
                weight=stan_W)
      if (type=="outer") {
        data$lambda_diag=as.array(lambda_diag)
        model=csnorm:::stanmodels$gauss_decay_outer
      } else {
        model=csnorm:::stanmodels$gauss_decay_perf
      }
      #optimize from scratch, to avoid getting stuck. Slower but more robust
      ops=optimize_stan_model(model=model, data=data, iter=cs@settings$iter,
                             verbose=verbose, init=0, init_alpha=init_alpha)
      
    } else {
      S_m2 = Diagonal(x=1/sdl^2)
      
      Xt = cbind(U_e,X)
      Dt = cbind(matrix(0,ncol=Dsets, nrow=Kdiag-2),D)
      tmp_X_S_m2_X = crossprod(Diagonal(x=1/sdl)%*%Xt)
      tmp_X_S_m2_k = t(Xt)%*%S_m2%*%kappa_hat
      DtD = crossprod(Dt)
      diags = list(rep(1,Kdiag), rep(-2,Kdiag))
      C=-bandSparse(Kdiag, Kdiag-1, k=c(0,-1),diagonals=list(diags[[1]],-diags[[1]]))
      Ct=rbind(matrix(0,nrow=Dsets,ncol=Kdiag), cbind(crossprod(X,W),C))
      
      epsilon = 1
      maxiter = 0
      
      while(epsilon > convergence_epsilon && maxiter < max_perf_iteration) {
        
        At = tmp_X_S_m2_X + Kdiag^2*lambda_diag^2*DtD
        
        At.PD = nearPD(At)$mat  
        fit = solve.QP(At.PD, tmp_X_S_m2_k, -Ct, meq = Dsets)
        betat = fit$solution
        
        eC=as.array(betat[1:Dsets])
        beta=betat[(Dsets+1):(Kdiag+Dsets)]
        
        if(type == 'outer') break
        
        nlambda_diag = (Kdiag - 2)/((Kdiag^2)*crossprod(D%*%beta)+1)
        nlambda_diag = sqrt(as.numeric(nlambda_diag))
        
        epsilon = abs(lambda_diag-nlambda_diag)
        lambda_diag = nlambda_diag
        maxiter = maxiter+1
   
      }
      
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
    }
    
    if(fit_model=='stan') {
      all_log_mean_counts = c(all_log_mean_counts,as.array(ops$par$log_mean_counts))
      all_log_decay = c(all_log_decay,ops$par$log_decay)
      all_beta_diag_centered = c(all_beta_diag_centered,ops$par$beta_diag_centered)
      d = 1
      valid_beta_diag = guarantee_beta_diag_increasing(ops$par$beta_diag)
      for (d2 in 1:TotalDsets) {
        if(XD[d2] == uXD) {
          all_beta_diag[d2,] = c(0,valid_beta_diag)
          all_beta_diag_diff[d2,] = ops$par$beta_diag_diff[d,]
          d = d + 1
        }
      }
      all_eC = c(all_eC,ops$par$eC)
      
      if (type=="perf") {
        all_lambda_diag = c(all_lambda_diag,ops$par$lambda_diag)
      }
      #update par slot
      #ops$par$value=ops$value
      all_value = ops$value
      
      ops$par$log_mean_counts=NULL
      
    } else {
      for (d2 in 1:TotalDsets) {
        if(XD[d2] == uXD) {
          all_beta_diag[d2,] = as.array(beta)
          all_beta_diag_diff[d2,] = as.array(D%*%beta)
          all_beta_diag_centered = c(all_beta_diag_centered,as.array(beta))
        }
      }
      all_log_mean_counts = c(all_log_mean_counts,as.array(log_mean_counts))
      all_log_decay = c(all_log_decay,log_decay)
      all_eC = c(all_eC,eC)
      if (type=="perf") {
        all_lambda_diag = c(all_lambda_diag,lambda_diag)
      }
      
      mus = sdl
      all_value = sum(dnorm(kappa_hat, mean = as.array(all_log_mean_counts), sd = as.array(mus), log = TRUE))
      
    }
    
  }
  op = list()
  attr(op,'par')
  for(nattr in list('beta_diag_centered','beta_diag','beta_diag_diff','log_decay','log_mean_counts','eC','lambda_diag','value','decay')) {
    attr(op$par,nattr)
  }
  op$par$beta_diag_centered=as.array(all_beta_diag_centered)
  op$par$beta_diag=as.matrix(all_beta_diag)
  op$par$beta_diag_diff=as.matrix(all_beta_diag_diff)
  op$par$log_mean_counts=as.array(all_log_mean_counts)
  op$par$log_decay=as.array(all_log_decay)
  op$par$eC=as.array(all_eC)
  if (type=="perf") {
    op$par$lambda_diag = as.array(all_lambda_diag)
  }
  op$par$value = all_value
  #make decay data table, reused at next call
  dmat=csd[,.(name,dbin,distance,kappahat,std,ncounts=weight,kappa=all_log_mean_counts)]
  setkey(dmat,name,dbin)
  op$par$decay=dmat 
  #rewrite log_decay as if it were calculated for each count
  dbins=cs@settings$dbins
  csub=cs@counts[,.(name,id1,id2,dbin=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  csd[,log_decay:=all_log_decay]
  a=csd[csub,.(name,id1,id2,log_decay),on=key(csd)]
  setkeyv(a,key(cs@counts))
  op$par$log_decay=a[,log_decay]
  #update par slot
  cs@par=modifyList(cs@par, op$par)
  return(cs)
}

#' Compute genomic etahat and weights using data as initial guess
#' @keywords internal
#' 
csnorm_gauss_genomic_muhat_data = function(cs, zbias, pseudocount=1e-2) {
  dispersion=cs@par$alpha
  #biases
  bts=rbind(cs@biases[,.(name,id,pos,cat="dangling L", etahat=log(dangling.L+pseudocount), std=sqrt(1/(dangling.L+pseudocount)+1/dispersion))],
            cs@biases[,.(name,id,pos,cat="dangling R", etahat=log(dangling.R+pseudocount), std=sqrt(1/(dangling.R+pseudocount)+1/dispersion))],
            cs@biases[,.(name,id,pos,cat="rejoined", etahat=log(rejoined+pseudocount), std=sqrt(1/(rejoined+pseudocount)+1/dispersion))])
  setkey(bts,id,name,cat)
  stopifnot(bts[,.N]==3*cs@biases[,.N])
  #counts
  cts=rbind(cs@counts[contact.close>0,.(name,id=id1,pos=pos1, cat="contact R", count=contact.close, weight=1)],
            cs@counts[contact.far>0,  .(name,id=id1,pos=pos1, cat="contact L", count=contact.far, weight=1)],
            cs@counts[contact.down>0, .(name,id=id1,pos=pos1, cat="contact R", count=contact.down, weight=1)],
            cs@counts[contact.up>0,   .(name,id=id1,pos=pos1, cat="contact L", count=contact.up, weight=1)],
            cs@counts[contact.far>0,  .(name,id=id2,pos=pos2, cat="contact R", count=contact.far, weight=1)],
            cs@counts[contact.close>0,.(name,id=id2,pos=pos2, cat="contact L", count=contact.close, weight=1)],
            cs@counts[contact.down>0, .(name,id=id2,pos=pos2, cat="contact R", count=contact.down, weight=1)],
            cs@counts[contact.up>0,   .(name,id=id2,pos=pos2, cat="contact L", count=contact.up, weight=1)],
            zbias[,.(name, id, pos, cat, count=0, weight=nzero)])
  cts[,c("etahat","var"):=list(log(count+pseudocount),(1/(count+pseudocount)+1/dispersion))]
  cts=cts[,.(etahat=weighted.mean(etahat,weight/var),std=sqrt(2/sum(weight/var)),weight=sum(weight)),
          keyby=c("id","pos","name","cat")]
  cts[,weight:=NULL] #not needed. Weights differ slightly because of dmin cutoff but it doesnt matter much
  setkeyv(cts,c("id","name","cat"))
  stopifnot(cts[,.N]==2*cs@biases[,.N])
  return(list(bts=bts,cts=cts))
}

#' Compute genomic etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_genomic_muhat_mean = function(cs, zbias) {
  #compute bias means
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
  bsub=bsub[,.(id,log_iota,log_rho)]
  #add bias informations to positive counts
  csub=copy(cs@counts)
  csub[,c("distance","log_decay"):=list(NULL,init$log_decay)]
  csub=merge(bsub[,.(id1=id,log_iota,log_rho)],csub,by="id1",all.x=F,all.y=T)
  csub=merge(bsub[,.(id2=id,log_iota,log_rho)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  csub=merge(cbind(cs@design[,.(name)],eC=init$eC), csub, by="name",all.x=F,all.y=T)
  #add it to zero counts
  zeta = merge(bsub, zbias, by="id")
  zeta = merge(cbind(cs@design[,.(name)],eC=init$eC), zeta, by="name",all.x=F,all.y=T)
  zeta[cat=="contact L",eta:=eC+log_iota]
  zeta[cat=="contact R",eta:=eC+log_rho]
  rm(bsub)
  #compute means
  csub[,lmu.base:=eC + log_decay]
  csub[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(lmu.base+log_iota1+log_rho2,
                                                           lmu.base+log_rho1 +log_rho2,
                                                           lmu.base+log_rho1 +log_iota2,
                                                           lmu.base+log_iota1+log_iota2)]
  csub[,lmu.base:=NULL]
  #collect all counts on left/right side
  cts=rbind(csub[contact.close>0,.(name, id=id1, pos=pos1, cat="contact R", count=contact.close, mu=exp(lmu.close), eta=eC + log_rho1,  weight=1)],
            csub[contact.far>0,  .(name, id=id1, pos=pos1, cat="contact L", count=contact.far,   mu=exp(lmu.far),   eta=eC + log_iota1, weight=1)],
            csub[contact.down>0, .(name, id=id1, pos=pos1, cat="contact R", count=contact.down,  mu=exp(lmu.down),  eta=eC + log_rho1,  weight=1)],
            csub[contact.up>0,   .(name, id=id1, pos=pos1, cat="contact L", count=contact.up,    mu=exp(lmu.up),    eta=eC + log_iota1, weight=1)],
            csub[contact.far>0,  .(name, id=id2, pos=pos2, cat="contact R", count=contact.far,   mu=exp(lmu.far),   eta=eC + log_rho2,  weight=1)],
            csub[contact.close>0,.(name, id=id2, pos=pos2, cat="contact L", count=contact.close, mu=exp(lmu.close), eta=eC + log_iota2, weight=1)],
            csub[contact.down>0, .(name, id=id2, pos=pos2, cat="contact R", count=contact.down,  mu=exp(lmu.down),  eta=eC + log_rho2,  weight=1)],
            csub[contact.up>0,   .(name, id=id2, pos=pos2, cat="contact L", count=contact.up,    mu=exp(lmu.up),    eta=eC + log_iota2, weight=1)],
            zeta[,.(name,id,pos,cat, count=0, mu=exp(eta), eta, weight=nzero)])
  cts[,var:=1/mu+1/init$alpha]
  cts=cts[,.(etahat=weighted.mean(count/mu-1+eta, weight/var),std=sqrt(2/sum(weight/var))), by=c("name","id","pos","cat")]
  setkey(cts,id,name,cat)
  stopifnot(cts[,.N]==2*cs@biases[,.N])
  return(list(bts=bts,cts=cts))
}


#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
csnorm_gauss_genomic = function(cs, zbias, verbose=T, init.mean="mean", init_alpha=1e-7, type=c("perf","outer"),fit_model=c('stan','nostan'), max_perf_iteration=1000, convergence_epsilon=1e-5) {
  type=match.arg(type)
  fit_model=match.arg(fit_model)
  if (init.mean=="mean") {
    a = csnorm:::csnorm_gauss_genomic_muhat_mean(cs, zbias)
  } else {
    a = csnorm_gauss_genomic_muhat_data(cs, zbias)
  }
  bts=a$bts
  cts=a$cts
  XB = as.array(cs@design[,genomic])
  
  #run optimization
  Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
  Totalbbegin=c(1,cs@biases[,.(name,row=.I)][name!=shift(name),row],cs@biases[,.N+1])
  TotalDsets = cs@design[,.N]
  
  all_beta_iota = c()
  all_beta_rho = c()
  all_beta_iota_diff = matrix(0, nrow = TotalDsets, ncol = Krow-2)
  all_beta_rho_diff = matrix(0, nrow = TotalDsets, ncol = Krow-2)
  all_log_iota = c()
  all_log_rho = c()
  all_log_mean_RJ = c()
  all_log_mean_DL = c()
  all_log_mean_DR = c()
  all_log_mean_cleft  = c()
  all_log_mean_cright = c()
  all_eC = c()
  all_eDE = c()
  all_eRJ = c()
  if (type=="perf") {
    all_lambda_iota = c()
    all_lambda_rho = c()
  }
  
  for(uXB in unique(XB)) {
    Dsets = 0
    bbegin = c()
    bbegin_stan = c(1)
    stanXB = c()
    for (d in 1:TotalDsets) {
      if(XB[d] == uXB) {
        bbegin = c(bbegin,Totalbbegin[d],Totalbbegin[d+1])
        bbegin_stan = c(bbegin_stan,(Totalbbegin[d+1]-Totalbbegin[d]+tail(bbegin_stan, n=1)))
        Dsets = Dsets + 1
        stanXB = c(stanXB,1) 
      }
    }
    SD = bbegin[2]-bbegin[1]
    #sparse spline design
    splinedegree=3 #Cubic spline
    cutsites = cs@biases[,pos][bbegin[1]:(bbegin[2]-1)]
    dx = 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree)
    t = min(cutsites) - dx*0.01 + dx * seq(-splinedegree,Krow-splinedegree+3)
    Bsp = spline.des(cutsites, knots = t, outer.ok = T, sparse=T)$design
    X = rbind(cbind(Bsp/2,Bsp/2),bdiag(Bsp,Bsp),bdiag(Bsp,Bsp))
    W=Matrix(cbind(c(rep(0,3*SD),rep(1,SD),rep(0,SD)), c(rep(0,4*SD),rep(1,SD))), sparse = T)
    diags = list(rep(1,Krow), rep(-2,Krow))
    D1 = bandSparse(Krow-2, Krow, k=0:2, diagonals=c(diags, diags[1]))
    if (type=="outer" || (!is.null(cs@par$lambda_iota) && !is.null(cs@par$lambda_rho))) {
      lambda_iota = cs@par$lambda_iota[uXB]
      lambda_rho = cs@par$lambda_rho[uXB]
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
        ncutsites = cs@biases[,pos][bbegin[(2*d-1)]:(bbegin[2*d]-1)]
        dx = 1.01*(max(ncutsites)-min(ncutsites))/(Krow-splinedegree)
        t = min(ncutsites) - dx*0.01 + dx * seq(-splinedegree,Krow-splinedegree+3)
        nBsp = spline.des(ncutsites, knots = t, outer.ok = T, sparse=T)$design
        cutsites = c(cutsites,ncutsites)
        nX = rbind(cbind(nBsp/2,nBsp/2),bdiag(nBsp,nBsp),bdiag(nBsp,nBsp))
        SDd = bbegin[2*d]-bbegin[(2*d-1)]
        X = rbind(X,nX)
        Bsp = bdiag(Bsp,nBsp)
        W = rbind(W,cbind(rep(0,5*SDd),rep(0,5*SDd)))
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
    if(fit_model=='stan') {
      
      data=list(DDsets=Dsets, Biases=1, XB=as.array(stanXB),
                Krow=Krow, SD=SD, bbegin=as.array(bbegin_stan), cutsitesD=cutsites,
                eta_hat_RJ=eta_hat_RJ, sd_RJ=sd_RJ,
                eta_hat_DL=eta_hat_DL, sd_DL=sd_DL,
                eta_hat_DR=eta_hat_DR, sd_DR=sd_DR,
                eta_hat_L=eta_hat_L, sd_L=sd_L,
                eta_hat_R=eta_hat_R, sd_R=sd_R)
      
      if (type=="outer") {
        model=csnorm:::stanmodels$gauss_genomic_outer
        data$lambda_rho=as.array(lambda_rho)
        data$lambda_iota=as.array(lambda_iota)
      } else {
        model=csnorm:::stanmodels$gauss_genomic_perf
      }
      
      op=optimize_stan_model(model=model, data=data, iter=cs@settings$iter,
                             verbose=verbose, init=0, init_alpha=init_alpha)
      
    } else {
    
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
      
      epsilon = 1
      maxiter = 0
      D = bdiag(lambda_iota*D1,lambda_rho*D1)
      DtD = crossprod(D)
      cholA = Cholesky(tmp_X_S_m2_X + Krow^2*DtD)
      while(epsilon > convergence_epsilon && maxiter < max_perf_iteration) {
        tmp_WtXAm1 = t(solve(cholA,tmp_Xt_W)) #2xK
        
        Gamma_v = solve(tmp_WtXAm1 %*% tmp_Xt_W, tmp_WtXAm1) #2xK
        #all(abs((Diagonal(2*Krow)-tmp_Xt_W%*%Gamma_v)%*%(Diagonal(2*Krow)-tmp_Xt_W%*%Gamma_v)-(Diagonal(2*Krow)-tmp_Xt_W%*%Gamma_v))<1e-5)
        
        tmp_Xt_Sm2_etas = crossprod(X,S_m2%*%etas) #Kx1
        beta_y = solve(cholA, tmp_Xt_Sm2_etas)
        beta_y = beta_y - solve(cholA, tmp_Xt_W) %*% (Gamma_v %*% tmp_Xt_Sm2_etas)
        
        tmp_Xt_Sm2_U = crossprod(X,S_m2%*%U_e)
        beta_U = solve(cholA, tmp_Xt_Sm2_U)
        beta_U = beta_U - solve(cholA, tmp_Xt_W) %*% (Gamma_v %*% tmp_Xt_Sm2_U)
        
        e=solve(t(U_e)%*%S_m2%*%(U_e-X%*%beta_U),t(U_e)%*%S_m2%*%(etas-X%*%beta_y))
        
        beta = beta_y - beta_U%*%e
        beta_iota = beta[1:Krow]
        beta_rho = beta[(Krow+1):(2*Krow)]
        
        if(type == 'outer') break
        
        nlambda_iota = (Krow - 2)/((Krow**2)*crossprod(D1%*%beta_iota)+1e6)
        nlambda_iota = sqrt(as.numeric(nlambda_iota))
        nlambda_rho = (Krow - 2)/((Krow**2)*crossprod(D1%*%beta_rho)+1e6)
        nlambda_rho = sqrt(as.numeric(nlambda_rho))
        
        epsilon = max(abs(lambda_iota-nlambda_iota),abs(lambda_rho-nlambda_rho))
        lambda_iota = nlambda_iota
        lambda_rho = nlambda_rho
        maxiter = maxiter+1
        
        D = bdiag(lambda_iota*D1,lambda_rho*D1)
        DtD = crossprod(D)
        cholA = update(cholA,tmp_X_S_m2_X + Krow^2*DtD)

      }
      
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
    }
    
    if(fit_model=='stan') {
      all_beta_iota = c(all_beta_iota,as.array(op$par$beta_iota))
      all_beta_rho = c(all_beta_rho,as.array(op$par$beta_rho))
      d = 1
      for (d2 in 1:TotalDsets) {
        if(XB[d2] == uXB) {
          all_beta_iota_diff[d2,] = op$par$beta_iota_diff[d,]
          all_beta_rho_diff[d2,] = op$par$beta_rho_diff[d,]
          d = d + 1
        }
      }
      all_log_iota = c(all_log_iota,as.array(op$par$log_iota))
      all_log_rho = c(all_log_rho,as.array(op$par$log_rho))
      all_log_mean_RJ = c(all_log_mean_RJ,as.array(op$par$log_mean_RJ))
      all_log_mean_DL = c(all_log_mean_DL,as.array(op$par$log_mean_DL))
      all_log_mean_DR = c(all_log_mean_DR,as.array(op$par$log_mean_DR))
      all_log_mean_cleft  = c(all_log_mean_cleft,as.array(op$par$log_mean_cleft))
      all_log_mean_cright = c(all_log_mean_cright,as.array(op$par$log_mean_cright))
      all_eC = c(all_eC,op$par$eC)
      all_eRJ = c(all_eRJ,op$par$eRJ)
      all_eDE = c(all_eDE,op$par$eDE)
      if (type=="perf") {
        all_lambda_iota = c(all_lambda_iota,op$par$lambda_iota)
        all_lambda_rho = c(all_lambda_rho,op$par$lambda_rho)
      }
      #update par slot
      op$par$value=op$value
      op$par$log_mean_DL=NULL
      op$par$log_mean_DR=NULL
      op$par$log_mean_RJ=NULL
      op$par$log_mean_cleft=NULL
      op$par$log_mean_cright=NULL
      
    } else {
      # make beta compatible with stan dispersion next step
      # beta_iota = beta_iota - beta_iota[1]
      # beta_iota = beta_iota[2:Krow]
      # beta_rho = beta_rho - beta_rho[1]
      # beta_rho = beta_rho[2:Krow]
      all_beta_iota = c(all_beta_iota,as.array(beta_iota))
      all_beta_rho = c(all_beta_rho,as.array(beta_rho))
      for (d2 in 1:TotalDsets) {
        if(XB[d2] == uXB) {
          all_beta_iota_diff[d2,] = as.array(D1%*%beta_iota)
          all_beta_rho_diff[d2,] = as.array(D1%*%beta_rho)
        }
      }
      all_log_iota = c(all_log_iota,as.array(log_iota))
      all_log_rho = c(all_log_rho,as.array(log_rho))
      all_log_mean_RJ = c(all_log_mean_RJ,as.array(log_mean_RJ))
      all_log_mean_DL = c(all_log_mean_DL,as.array(log_mean_DL))
      all_log_mean_DR = c(all_log_mean_DR,as.array(log_mean_DR))
      all_log_mean_cleft  = c(all_log_mean_cleft,as.array(log_mean_cleft))
      all_log_mean_cright = c(all_log_mean_cright,as.array(log_mean_cright))
      all_eC = c(all_eC,eC)
      all_eRJ = c(all_eRJ,eRJ)
      all_eDE = c(all_eDE,eDE)
      if (type=="perf") {
        all_lambda_iota = c(all_lambda_iota,lambda_iota)
        all_lambda_rho = c(all_lambda_rho,lambda_rho)
      }
      op = list(par=list())
      op$par$beta_iota=as.array(all_beta_iota)
      op$par$beta_rho=as.array(all_beta_rho)
      op$par$beta_rho_diff=as.matrix(all_beta_rho_diff)
      op$par$beta_iota_diff=as.matrix(all_beta_iota_diff)
      op$par$log_iota=as.array(all_log_iota)
      op$par$log_rho=as.array(all_log_rho)
      op$par$eC=as.array(all_eC)
      op$par$eRJ=as.array(all_eRJ)
      op$par$eDE=as.array(all_eDE)
      if (type=="perf") {
        op$par$lambda_iota = as.array(all_lambda_iota)
        op$par$lambda_rho = as.array(all_lambda_rho)
      }
      means = cbind(all_log_mean_RJ,all_log_mean_DL,all_log_mean_DR,all_log_mean_cleft,all_log_mean_cright)
      mus = cbind(sd_RJ,sd_DL,sd_DR,sd_L,sd_R)
      op$par$value = sum(dnorm(etas, mean = as.array(means), sd = as.array(mus), log = TRUE))
    }
  }
  
  #make nice output table
  bout=rbind(bts[cat=="dangling L",.(cat, name, id, pos, etahat, std, eta=as.array(all_log_mean_DL))],
             bts[cat=="dangling R",.(cat, name, id, pos, etahat, std, eta=as.array(all_log_mean_DR))],
             bts[cat=="rejoined",.(cat, name, id, pos, etahat, std, eta=as.array(all_log_mean_RJ))])
  cout=rbind(cts[cat=="contact L",.(cat, name, id, pos, etahat, std, eta=as.array(all_log_mean_cleft))],
             cts[cat=="contact R",.(cat, name, id, pos, etahat, std, eta=as.array(all_log_mean_cright))])
  bout=rbind(bout,cout)
  setkey(bout, id, name, cat)
  op$par$biases=bout
  cs@par=modifyList(cs@par, op$par)
  gc()
  return(cs)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_gauss_dispersion = function(cs, counts, weight=cs@design[,.(name,wt=1)], verbose=T, init_alpha=1e-7, type=c("outer","perf"), ncores=ncores) {
  type=match.arg(type)
  #predict all means and put into table
  counts=csnorm:::csnorm_predict_all_parallel(cs,counts,ncores = ncores)
  setkeyv(counts,key(cs@counts))
  stopifnot(cs@biases[,.N]==length(cs@par$log_iota))
  #
  #fit dispersion and exposures
  bbegin=c(1,cs@biases[,.(name,row=.I)][name!=shift(name),row],cs@biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=cs@design[,.N], Biases=cs@design[,uniqueN(genomic)], Decays=cs@design[,uniqueN(decay)],
               XB=as.array(cs@design[,genomic]), XD=as.array(cs@design[,decay]),
               SD=cs@biases[,.N], bbegin=bbegin,
               cutsitesD=cs@biases[,pos], rejoined=cs@biases[,rejoined],
               danglingL=cs@biases[,dangling.L], danglingR=cs@biases[,dangling.R],
               N=counts[,.N], cbegin=cbegin,
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight[,wt]),
               log_iota=cs@par$log_iota, log_rho=cs@par$log_rho,
               log_mean_cclose=counts[,log_mean_cclose], log_mean_cfar=counts[,log_mean_cfar],
               log_mean_cup=counts[,log_mean_cup], log_mean_cdown=counts[,log_mean_cdown])
  init=list(eC_sup=as.array(counts[,log(mean(contact.close/exp(log_mean_cclose))),by=name][,V1]),
            eRJ=as.array(cs@biases[,.(name,frac=rejoined/exp((cs@par$log_iota+cs@par$log_rho)/2))][,log(mean(frac)),by=name][,V1]),
            eDE=as.array(cs@biases[,.(name,frac=(dangling.L/exp(cs@par$log_iota)+dangling.R/exp(cs@par$log_rho))/2)][
              ,log(mean(frac)),by=name][,V1]))
  init$mu=mean(exp(init$eC_sup[1]+counts[name==name[1],log_mean_cclose]))
  init$alpha=max(0.001,1/(var(counts[name==name[1],contact.close]/init$mu)-1/init$mu))
  init$mu=NULL
  op=optimize_stan_model(model=csnorm:::stanmodels$gauss_dispersion, data=data, iter=cs@settings$iter,
                         verbose=verbose, init=init, init_alpha=init_alpha)
  cs@par=modifyList(cs@par, op$par[c("eRJ","eDE","alpha")])
  cs@par$eC=cs@par$eC+op$par$eC_sup
  #
  #estimate lambdas if outer iteration
  Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
  Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))* cs@settings$bf_per_decade)
  if (type=="outer") {
    lambdas = copy(cs@design)
    lambdas[,lambda_iota:=sqrt((Krow-2)/(Krow^2*diag(tcrossprod(cs@par$beta_iota_diff))+1e6))] #sigma=1e-3 for genomic
    lambdas[,lambda_rho:=sqrt((Krow-2)/(Krow^2*diag(tcrossprod(cs@par$beta_rho_diff))+1e6))]
    lambdas[,lambda_diag:=sqrt((Kdiag-2)/(Kdiag^2*diag(tcrossprod(cs@par$beta_diag_diff))+1))] #sigma=1 for decay
    
    cs@par$lambda_iota=as.array(lambdas[unique(genomic)][order(genomic),lambda_iota])
    cs@par$lambda_rho=as.array(lambdas[unique(genomic)][order(genomic),lambda_rho])
    cs@par$lambda_diag=as.array(lambdas[unique(decay)][order(decay),lambda_diag])
  }
  #
  #compute log-posterior
  cs@par$value = op$value + (Krow-2)/2*sum(log(cs@par$lambda_iota/exp(1))+log(cs@par$lambda_rho/exp(1))) +
                            (Kdiag-2)/2*sum(log(cs@par$lambda_diag/exp(1)))
  return(cs)
}


#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
has_converged = function(cs) {
  params=cs@diagnostics$params
  if (params[,.N]<=3) return(FALSE)
  tol=cs@settings$tol.obj
  laststep=params[,step[.N]]
  delta=abs(params[step==laststep,value]-params[step==laststep-1,value])
  return(all(delta<tol))
}

#' count number of zeros in each decay bin
#' @keywords internal
#' 
get_nzeros_per_decay = function(cs, ncores=1) {
  #bin distances
  dbins=cs@settings$dbins
  stopifnot(cs@counts[id1>=id2,.N]==0)
  mcounts=melt(cs@counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0]
  mcounts[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #Count positive counts in these bins
  Nkd = mcounts[,.(nnz=.N),keyby=c("name","bdist")]
  #Nkd[,mdist:=sqrt(dbins[unclass(bdist)+1]*dbins[unclass(bdist)])] #cannot use dbins because some rows could be missing
  #stopifnot(Nkd[mdist<cs@settings$dmin,.N]==0) #otherwise cs@counts has not been censored properly
  #Count the number of crossings per distance bin
  #looping over IDs avoids building NxN matrix
  registerDoParallel(cores=ncores)
  Nkz = foreach(i=cs@biases[,(1:.N)], .multicombine = T,
                .combine=function(...){rbind(...)[,.(ncross=sum(ncross)),keyby=c("name","bdist")]}) %dopar% {
    stuff=c(cs@biases[i,.(name,id,pos)])
    dists=cs@biases[name==stuff$name & id>stuff$id,.(name,distance=abs(pos-stuff$pos))]
    if (cs@settings$circularize>0)  dists[,distance:=pmin(distance,cs@settings$circularize+1-distance)]
    dists=dists[distance>=cs@settings$dmin]
    dists[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
    dists[,.(ncross=.N),keyby=c("name","bdist")]
  }
  #deduce zero counts
  Nkz = merge(Nkz,Nkd,all=T)
  Nkz[is.na(nnz),nnz:=0]
  Nkz[,nzero:=4*ncross-nnz]
  return(Nkz)
}

#' count number of zeros in each decay bin
#' @keywords internal
#' 
get_nzeros_per_cutsite = function(cs, ncores=1) {
  stopifnot(cs@counts[id1>=id2,.N]==0)
  cts=rbind(cs@counts[contact.close>0,.(name, id=id1, cat="contact R", count=contact.close)],
            cs@counts[contact.far>0,  .(name, id=id1, cat="contact L", count=contact.far)],
            cs@counts[contact.down>0, .(name, id=id1, cat="contact R", count=contact.down)],
            cs@counts[contact.up>0,   .(name, id=id1, cat="contact L", count=contact.up)],
            cs@counts[contact.far>0,  .(name, id=id2, cat="contact R", count=contact.far)],
            cs@counts[contact.close>0,.(name, id=id2, cat="contact L", count=contact.close)],
            cs@counts[contact.down>0, .(name, id=id2, cat="contact R", count=contact.down)],
            cs@counts[contact.up>0,   .(name, id=id2, cat="contact L", count=contact.up)])
  cts=cts[,.N,keyby=c("id","name","cat")]
  zbias=rbind(cts[cat=="contact L"][cs@biases,.(name,id,pos,cat="contact L",nnz=N)],
              cts[cat=="contact R"][cs@biases,.(name,id,pos,cat="contact R",nnz=N)])
  zbias[is.na(nnz),nnz:=0]
  setkey(zbias,id,name)
  #count masked contacts where d<dmin
  registerDoParallel(cores=ncores)
  masked = foreach(i=cs@biases[,(1:.N)], .combine=rbind) %dopar% {
    stuff=c(cs@biases[i,.(name,id,pos)])
    dists=cs@biases[name==stuff$name & id!=stuff$id,.(name,distance=abs(pos-stuff$pos))]
    if (cs@settings$circularize>0)  dists[,distance:=pmin(distance,cs@settings$circularize+1-distance)]
    dists[distance<cs@settings$dmin,.(name=name[1],id=stuff$id,ncross.close=.N)]
  }
  setkey(masked,id,name)
  masked=masked[cs@biases,.(name,id,ncross.close)]
  masked[is.na(ncross.close),ncross.close:=0]
  #deduce number of zero counts
  zbias=masked[zbias]
  zbias[,ncs:=(.N/2-1),by=name]
  zbias[,nzero:=2*(ncs-ncross.close)-nnz]
  stopifnot(zbias[,.N]==cs@biases[,.N*2])
  setkey(zbias,id,name,cat)
  return(zbias)
}

#' Get the first ncounts/d of each of d datasets, including zeros 
#' 
#' count is always a bit smaller because we censor those that are <dmin without adding more counts
#' 
#' @keywords internal
#' 
subsample_counts = function(cs, ncounts, dset=NA) {
  ncounts_per_dset=as.integer(ncounts/cs@design[,.N])
  if (is.na(dset)) {
    cts = foreach (d=cs@design[,name],.combine=rbind) %do% subsample_counts(cs,ncounts_per_dset,dset=d)
  } else {
    #get name and id of counts
    nbiases = cs@biases[name==dset, .N]
    ncounts = min(nbiases*(nbiases-1)/2,ncounts)
    ids=c(cs@biases[name==dset,.(minid=min(id),maxid=max(id))])
    cts=data.table(name=dset,id1=ids$minid,id2=(ids$minid+1):ids$maxid)
    while(cts[,.N]<ncounts)
      cts=rbind(cts,data.table(name=dset,id1=cts[.N,id1+1],id2=cts[.N,id1+2]:ids$maxid))
    cts=cts[1:ncounts]
    #merge positions and compute distances
    cts = merge(cts,cs@biases[,.(name,id,pos)],by.x=c("name","id1"),by.y=c("name","id"))
    cts = merge(cts,cs@biases[,.(name,id,pos)],by.x=c("name","id2"),by.y=c("name","id"),suffixes=c("1","2"))
    cts[,distance:=abs(pos2-pos1)]
    if (cs@settings$circularize>0) cts[,distance:=pmin(distance, cs@settings$circularize+1-distance)] 
    cts = cts[distance>=cs@settings$dmin]
    #merge counts
    setkey(cts, id1, id2, name)
    cts = merge(cts, cs@counts[,.(name,id1,id2,contact.close,contact.down,contact.far,contact.up)], all.x=T)
    cts[is.na(contact.close),contact.close:=0]
    cts[is.na(contact.down),contact.down:=0]
    cts[is.na(contact.far),contact.far:=0]
    cts[is.na(contact.up),contact.up:=0]
    if (cts[,uniqueN(c(contact.close,contact.far,contact.up,contact.down))]<2)
      stop("dataset too sparse, please increase ncounts")
    return(cts)
  }
}

#' Diagnostics plots to monitor convergence of normalization (gaussian
#' approximation)
#' 
#' @param cs a normalized CSnorm object
#'   
#' @return Two plots in a list. The first is for the three log-likelihoods, the
#'   second is for the parameters.
#' @export
#' 
#' @examples
plot_diagnostics = function(cs) {
  plot=ggplot(cs@diagnostics$params[,.(step,leg,value,out.last)])+
    geom_line(aes(step,value))+geom_point(aes(step,value,colour=out.last))+facet_wrap(~leg, scales = "free")+
    theme(legend.position="bottom")
  vars=foreach(var=c("eC","eRJ","eDE","alpha","lambda_iota","lambda_rho","lambda_diag"),
               trans=(c("exp","exp","exp",NA,"log","log","log")),.combine=rbind) %do% get_all_values(cs,var,trans)
  plot2=ggplot(vars)+geom_line(aes(step,value))+
    geom_point(aes(step,value,colour=leg))+facet_wrap(~variable, scales = "free_y")
  return(list(plot=plot,plot2=plot2))
}

#' Cut-site normalization (normal approximation)
#' 
#' Alternates two approximations to the exact model, fitting the diagonal decay
#' and iota/rho.
#' 
#' @param cs CSnorm object as returned by \code{\link{merge_cs_norm_datasets}}
#' @param bf_per_kb positive numeric. Number of cubic spline basis functions per
#'   kilobase, for genomic bias estimates. Small values make the optimization 
#'   easy, but makes the genomic biases stiffer.
#' @param bf_per_decade positive numeric. Number of cubic spline basis functions
#'   per distance decade (in bases), for diagonal decay. Default parameter 
#'   should suffice.
#' @param bins_per_bf positive integer. Number of distance bins to split basis 
#'   functions into. Must be sufficiently small so that the diagonal decay is 
#'   approximately constant in that bin.
#' @param lambdas positive numeric. Length scales to try out as initial
#'   condition.
#' @param ngibbs positive integer. Number of gibbs sampling iterations.
#' @param iter positive integer. Number of optimization steps for each stan 
#'   optimization call.
#' @param fit.decay,fit.genomic,fit.disp boolean. Whether to fit diagonal decay or
#'   genomic biases. Set to FALSE only for diagnostics.
#' @param verbose Display progress if TRUE
#' @param init_alpha positive numeric, default 1e-5. Initial step size of LBFGS
#'   line search.
#' @param init.dispersion positive numeric. Value of the dispersion to use initially.
#' @param tol.obj positive numeric (default 1e-1). Convergence tolerance on changes in the three likelihoods.
#' @param type character. Type of iteration to perform. The outer iteration ("outer", default) optimizes
#'   penalties with the dispersion parameter. The performance iteration ("perf") optimizes penalties with
#'   the biases. It is recommended to try performance iteration, but if convergence fails to revert to outer
#'   iteration.
#' @param ncores positive integer (default 1). Number of cores to use for initialization and dispersion step.
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
#' 
run_gauss = function(cs, init=NULL, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10,
                     ngibbs = 3, iter=10000, fit.decay=T, fit.genomic=T, fit.disp=T,
                     verbose=T, ncounts=100000, init_alpha=1e-7, init.dispersion=10,
                     tol.obj=1e-1, type="outer",fit_model="nostan", ncores=1) {
  type=match.arg(type,c("outer","perf"))
  fit_model=match.arg(fit_model,c("stan","nostan"))
  if (verbose==T) cat("Normalization with fast approximation and ",type,"iteration\n")
  #clean object if dirty
  cs@par=list() #in case we have a weird object
  cs@binned=list()
  setkey(cs@biases, id, name)
  setkey(cs@counts, id1, id2, name)
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #add settings
  cs@settings = c(cs@settings[c("circularize","dmin","dmax","qmin","qmax")],
                  list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
                       iter=iter, init_alpha=init_alpha, init.dispersion=init.dispersion, tol.obj=tol.obj))
  #get number of zeros along cut sites and decay
  stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
  cs@settings$dbins=10**seq(log10(cs@settings$dmin),log10(cs@settings$dmax)+stepsz,stepsz)
  if(verbose==T) cat("Counting zeros for decay\n")
  zdecay = csnorm:::get_nzeros_per_decay(cs, ncores=ncores)
  if(verbose==T) cat("Counting zeros for bias\n")
  zbias = csnorm:::get_nzeros_per_cutsite(cs, ncores=ncores)
  #
  if(verbose==T) cat("Subsampling counts for dispersion\n")
  subcounts = csnorm:::subsample_counts(cs, ncounts)
  subcounts.weight = merge(zdecay[,.(nc=sum(ncross)),by=name],subcounts[,.(ns=.N),keyby=name])[,.(name,wt=nc/ns)]
  #initial guess
  if (is.null(init)) {
    if (verbose==T) cat("No initial guess provided\n")
    cs@diagnostics=list()
    laststep=0
    init.mean="data"
    if (type=="outer") {
      cs=csnorm:::fill_parameters_outer(cs, dispersion=init.dispersion, fit.decay=fit.decay,
                               fit.genomic=fit.genomic, fit.disp=fit.disp)
    } else {
      cs=csnorm:::fill_parameters_perf(cs, dispersion=init.dispersion, fit.decay=fit.decay,
                               fit.genomic=fit.genomic, fit.disp=fit.disp)
    }
  } else {
    if (verbose==T) cat("Using provided initial guess\n")
    if (is.data.table(cs@diagnostics$params)) laststep = cs@diagnostics$params[,max(step)] else laststep = 0
    if (fit_model=="stan") init$beta_diag = guarantee_beta_diag_increasing(init$beta_diag)
    init.mean="mean"
    cs@par=init
  }
  #gibbs sampling
  for (i in (laststep + 1:ngibbs)) {
    #fit diagonal decay given iota and rho
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay ")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_decay(cs, zdecay, init.mean=init.mean,
                                                                               init_alpha=init_alpha, type=type,fit_model=fit_model)))
      if(length(output) == 0) { output = 'ok' }
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="decay", out=output, runtime=a[1]+a[4], type=type)
      if (verbose==T) cat("log-likelihood = ",cs@par$value, "\n")
    }
    #fit iota and rho given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic ")
      decay_eC = cs@par$eC
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_genomic(cs, zbias, init.mean=init.mean,
                                                                                 init_alpha=init_alpha, type=type,fit_model=fit_model)))
      if(length(output) == 0) { output = 'ok' }
      cs@par$eC = as.array(colMeans(rbind(decay_eC,cs@par$eC)))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="bias", out=output, runtime=a[1]+a[4], type=type)
      if (verbose==T) cat("log-likelihood = ",cs@par$value, "\n")
    }
    init.mean="mean"
    if (fit.disp==T) {
      #fit exposures and dispersion
      if (verbose==T) cat("Gibbs",i,": Remaining parameters ")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_dispersion(cs, counts=subcounts, weight=subcounts.weight,
                                                                                    init_alpha=init_alpha, type=type, ncores = ncores)))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="disp", out=output, runtime=a[1]+a[4], type=type)
      if (verbose==T) cat("log-likelihood = ",cs@par$value,"\n")
    }
    if (verbose==T) cat("Gibbs",i,": dispersion = ",cs@par$alpha, " lambda_iota = ",cs@par$lambda_iota, "\n")
    if (has_converged(cs)) break
  }
  if (verbose==T) cat("Done\n")
  if ("init" %in% names(init)) init$init=NULL
  cs@par$init=init
  return(cs)
}


