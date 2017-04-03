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
    Kdiag=cs@settings$Kdiag
    beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
    init=c(init,list(beta_diag=beta_diag, lambda_diag=array(lambda,dim=Decays)))
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
    Kdiag=cs@settings$Kdiag
    beta_diag=matrix(rep(seq(0.1,1,length.out = Kdiag-1), each=Decays), Decays, Kdiag-1)
    init=c(init,list(beta_diag=beta_diag))
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

#' Compute decay etahat and weights using data as initial guess. Signal is assumed to be zero.
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_data = function(cs, pseudocount=1e-2) {
  #bin distances
  dbins=cs@settings$dbins
  mcounts=melt(cs@counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0,.(name,distance,category,count)]
  mcounts[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #add zero counts
  zdecay = cs@zeros$bg[,.(nzero=sum(nzero)/2),keyby=c("name","dbin")] #each contact is counted twice
  mcounts = rbind(mcounts[,.(name,dbin,category,count,weight=1)],
                  zdecay[,.(name,dbin,category=NA,count=0,weight=nzero)])
  #compute z-scores and sum counts
  mcounts[,c("kappahat","var"):=list(log(count+pseudocount), 1/(count+pseudocount)+1/cs@par$alpha)]
  csd = mcounts[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
                   kappahat=weighted.mean(kappahat, weight/var),
                   std=1/sqrt(sum(weight/var)), weight=sum(weight)), keyby=c("name", "dbin")]
  stopifnot(csd[,!is.na(distance)])
  return(csd)
}

#' Compute means for positive and zero counts using previous params
#' @keywords internal
#' 
csnorm_gauss_common_muhat_mean = function(cs) {
  init=cs@par
  ### positive counts
  #compute means
  cpos=copy(cs@counts)
  bsub=cs@biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  cpos=merge(bsub[,.(id1=id,log_iota,log_rho)],cpos,by="id1",all.x=F,all.y=T)
  cpos=merge(bsub[,.(id2=id,log_iota,log_rho)],cpos,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  cpos=merge(cbind(cs@design[,.(name)],eC=init$eC), cpos, by="name",all.x=F,all.y=T)
  dbins=cs@settings$dbins
  cpos[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  cpos=merge(cpos,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  cpos[,log_mu.base:=eC + log_decay]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  cpos[,log_mu.base:=NULL]
  #collect all counts on left/right side, per distance
  cpos=rbind(cpos[contact.close>0,.(name, id=id1, pos=pos1, distance, cat="contact R",
                                    count=contact.close, lmu=lmu.close, weight=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.far>0,  .(name, id=id1, pos=pos1, distance, cat="contact L",
                                    count=contact.far,   lmu=lmu.far,   weight=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.down>0, .(name, id=id1, pos=pos1, distance, cat="contact R",
                                    count=contact.down,  lmu=lmu.down,  weight=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.up>0,   .(name, id=id1, pos=pos1, distance, cat="contact L",
                                    count=contact.up,    lmu=lmu.up,    weight=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.far>0,  .(name, id=id2, pos=pos2, distance, cat="contact R",
                                    count=contact.far,   lmu=lmu.far,   weight=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.close>0,.(name, id=id2, pos=pos2, distance, cat="contact L",
                                    count=contact.close, lmu=lmu.close, weight=1, eC, log_decay, log_bias=log_iota2)],
             cpos[contact.down>0, .(name, id=id2, pos=pos2, distance, cat="contact R",
                                    count=contact.down,  lmu=lmu.down,  weight=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.up>0,   .(name, id=id2, pos=pos2, distance, cat="contact L",
                                    count=contact.up,    lmu=lmu.up,    weight=1, eC, log_decay, log_bias=log_iota2)])
  dbins=cs@settings$dbins
  cpos[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  cpos[,distance:=NULL]
  ### zero counts
  czero = merge(init$decay[,.(name,dbin,log_decay)], cs@zeros$bg, by=c("name","dbin"))
  stopifnot(czero[is.na(log_decay),.N]==0)
  czero = merge(bsub,czero,by="id",all.x=F,all.y=T)
  czero[,log_bias:=ifelse(cat=="contact L",log_iota,log_rho)]
  czero = merge(cbind(cs@design[,.(name)],eC=init$eC), czero, by="name",all.x=F,all.y=T)
  czero[,lmu:=eC + log_decay + log_bias] #we don't average over j
  czero = czero[nzero>0,.(name,id,pos,dbin,cat,count=0,lmu,weight=nzero,log_bias,log_decay,eC)]
  cts=rbind(cpos,czero)
  ### add signal
  if (cs@zeros$sig[,.N]>0) {
    cts[,bin:=cut(pos, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    signal=merge(cs@zeros$sig, cs@par$signal, by=c("name","bin1","bin2"), all.x=T)
    stopifnot(all(complete.cases(signal)))
    signal=rbind(signal[,.(name,bin=bin1,dbin,ncross,phi)],signal[bin2>bin1,.(name,bin=bin2,dbin,ncross,phi)])
    signal=signal[,.(phi=weighted.mean(phi,ncross)),keyby=c("name","bin","dbin")]
    cts=signal[cts,,on=key(signal)]
    cts[,mu:=exp(lmu+phi)]
  } else {
    cts[,mu:=exp(lmu)]
  }
  ### finalize
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  #ggplot(cts[name=="T47D es 60 MboI 1"&cat=="contact R"])+geom_line(aes(dbin,log_decay,colour=count>0,group=count>0))
  cts[,c("count","mu","pos","lmu"):=NULL] #not needed downstream
  return(cts)
}

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_mean = function(cs) {
  cts = csnorm:::csnorm_gauss_common_muhat_mean(cs)
  cts[,kappaij:=eC+log_decay]
  dbins=cs@settings$dbins
  csd = cts[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
               kappahat=weighted.mean(z+kappaij, weight/var),
               std=1/sqrt(sum(weight/(2*var))), weight=sum(weight)/2), keyby=c("name", "dbin")] #each count appears twice
  return(csd)
}

#' Optimize decay parameters
#' @keywords internal
#' 
csnorm_gauss_decay_optimize = function(csd, design, Kdiag, original_lambda_diag, type, max_perf_iteration=1000, convergence_epsilon=1e-5) {
  Totalcbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  #cbegin=c(1,csd[,.(name,row=.I)][name!=shift(name),row],csd[,.N+1])
  TotalDsets = design[,.N]
  XD=as.array(design[,decay])
  
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
    for (d in 1:TotalDsets) {
      if(XD[d] == uXD) {
        cbegin = c(cbegin,Totalcbegin[d],Totalcbegin[d+1])
        Dsets = Dsets + 1
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
    if (type=="outer" || (!is.null(original_lambda_diag))) {
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
      #salt=Kdiag^2*lambda_diag^2*0.1*bdiag(Diagonal(Dsets)*0,Diagonal(Kdiag))
      #At.PD = At + salt
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
  dmat=csd[,.(name,dbin,distance,kappahat,std,ncounts=weight,kappa=all_log_mean_counts,log_decay=all_log_decay)]
  setkey(dmat,name,dbin)
  op$par$decay=dmat 
  return(op)
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
csnorm_gauss_decay = function(cs, init.mean="mean", type=c("outer","perf"), max_perf_iteration=1000, convergence_epsilon=1e-5) {
  type=match.arg(type)
  if (init.mean=="mean") {
    csd = csnorm:::csnorm_gauss_decay_muhat_mean(cs)
  } else {
    csd = csnorm:::csnorm_gauss_decay_muhat_data(cs)
  }
  #run optimization
  op = csnorm:::csnorm_gauss_decay_optimize(csd, cs@design, cs@settings$Kdiag, cs@par$lambda_diag,
                                            type, max_perf_iteration=max_perf_iteration, convergence_epsilon=convergence_epsilon)
  #update par slot
  cs@par=modifyList(cs@par, op$par)
  return(cs)
}

#' Compute genomic etahat and weights using data as initial guess
#' @keywords internal
#' 
csnorm_gauss_genomic_muhat_data = function(cs, pseudocount=1e-2) {
  dispersion=cs@par$alpha
  #biases
  bts=rbind(cs@biases[,.(name,id,pos,cat="dangling L", etahat=log(dangling.L+pseudocount), std=sqrt(1/(dangling.L+pseudocount)+1/dispersion))],
            cs@biases[,.(name,id,pos,cat="dangling R", etahat=log(dangling.R+pseudocount), std=sqrt(1/(dangling.R+pseudocount)+1/dispersion))],
            cs@biases[,.(name,id,pos,cat="rejoined", etahat=log(rejoined+pseudocount), std=sqrt(1/(rejoined+pseudocount)+1/dispersion))])
  setkey(bts,id,name,cat)
  stopifnot(bts[,.N]==3*cs@biases[,.N])
  #counts
  zbias = cs@zeros$bg[,.(nzero=sum(nzero)),keyby=c("id","name","cat")][cs@biases[,.(name,id,pos)]]
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

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_genomic_muhat_mean = function(cs) {
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
  cts.common = csnorm:::csnorm_gauss_common_muhat_mean(cs)
  cts.common[,etaij:=eC+log_bias]
  cts = cts.common[,.(etahat=weighted.mean(z+etaij, weight/var),std=sqrt(2/sum(weight/var))),
                   keyby=c("id","name","cat")]
  cts = merge(cts,cs@biases[,.(name,id,pos)],by=c("name","id"))
  setkey(cts, id, name, cat)
  return(list(bts=bts,cts=cts))
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
csnorm_gauss_genomic = function(cs, verbose=T, init.mean="mean", type=c("perf","outer"), max_perf_iteration=1000, convergence_epsilon=1e-5) {
  type=match.arg(type)
  if (init.mean=="mean") {
    a = csnorm:::csnorm_gauss_genomic_muhat_mean(cs)
  } else {
    a = csnorm:::csnorm_gauss_genomic_muhat_data(cs)
  }
  bts=a$bts
  cts=a$cts
  XB = as.array(cs@design[,genomic])
  
  #run optimization
  Krow=cs@settings$Krow
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
    for (d in 1:TotalDsets) {
      if(XB[d] == uXB) {
        bbegin = c(bbegin,Totalbbegin[d],Totalbbegin[d+1])
        Dsets = Dsets + 1
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
  #cs@par$eC=cs@par$eC+op$par$eC_sup
  #
  #estimate lambdas if outer iteration
  Krow=cs@settings$Krow
  Kdiag=cs@settings$Kdiag
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

#' fit signal using sparse fused lasso
#' @keywords internal
#' 
csnorm_gauss_signal = function(cs, verbose=T, ncores=ncores, tol=1e-3) {
  cts = csnorm:::csnorm_predict_binned_counts_irls(cs, cs@settings$base.res, cs@zeros$sig)
  mat = csnorm:::csnorm_compute_raw_signal(cts, cs@par$alpha, cs@par$signal)
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phihat))+facet_wrap(~name)
  #
  #perform fused lasso on signal
  groupnames=mat[,as.character(unique(name))]
  registerDoParallel(cores=ncores)
  params = foreach(g=groupnames, .combine=rbind) %dopar%
    csnorm:::csnorm_fused_lasso_signal(mat[name==g], cs@settings$trails, tol=tol, ncores=ncores, verbose=verbose)
  #compute matrix at new params
  #save(mat,params,file=paste0("mat_step_",step,".RData"))
  mat = foreach (g=groupnames, .combine=rbind) %dopar% {
    p=params[name==g]
    matg=mat[name==g]
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, cs@settings$trails, p$lambda1, p$lambda2, p$eCprime)]
    matg
  }
  #store new signal in cs and update eC
  mat[,phi:=value]
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phihat))+facet_wrap(~name)+scale_fill
  cs@par$signal=mat[,.(name,bin1,bin2,phi)]
  params=merge(cbind(cs@design[,.(name)],eC=cs@par$eC), params, by="name",all=T)
  cs@par$eC=as.array(params[,eC+eCprime])
  cs@par$eCprime=as.array(params[,eCprime])
  cs@par$lambda1=as.array(params[,lambda1])
  cs@par$lambda2=as.array(params[,lambda2])
  cs@par$value = params[,sum(BIC)]
  return(cs)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
has_converged = function(cs) {
  params=cs@diagnostics$params
  tol=cs@settings$tol.obj
  laststep=params[,step[.N]]
  delta=abs(params[step==laststep,value]-params[step==laststep-1,value])
  return(all(delta<tol))
}

#' count number of zeros in each decay and genomic bin
#' @keywords internal
#' 
get_nzeros_normalization = function(cs, ncores=1) {
  stopifnot(cs@counts[id1>=id2,.N]==0)
  #count left and right
  cts=rbind(cs@counts[contact.close>0,.(name, id=id1, distance, cat="contact R", count=contact.close)],
            cs@counts[contact.far>0,  .(name, id=id1, distance, cat="contact L", count=contact.far)],
            cs@counts[contact.down>0, .(name, id=id1, distance, cat="contact R", count=contact.down)],
            cs@counts[contact.up>0,   .(name, id=id1, distance, cat="contact L", count=contact.up)],
            cs@counts[contact.far>0,  .(name, id=id2, distance, cat="contact R", count=contact.far)],
            cs@counts[contact.close>0,.(name, id=id2, distance, cat="contact L", count=contact.close)],
            cs@counts[contact.down>0, .(name, id=id2, distance, cat="contact R", count=contact.down)],
            cs@counts[contact.up>0,   .(name, id=id2, distance, cat="contact L", count=contact.up)])
  #bin distances
  dbins=cs@settings$dbins
  cts[,dbin:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  #count per cut site and distance bin
  cts=cts[,.(nnz=.N),keyby=c("id","name","dbin","cat")]
  #Count the number of crossings per distance bin
  #looping over IDs avoids building NxN matrix
  registerDoParallel(cores=ncores)
  chunksize=cs@biases[,ceiling(.N/(10*ncores))]
  nchunks=cs@biases[,ceiling(.N/chunksize)]
  zeros = foreach(chunk=1:nchunks, .combine=rbind) %dopar% {
    bs=cs@biases[((chunk-1)*chunksize+1):min(.N,chunk*chunksize)]
    dists = foreach(n=bs[,name], i=bs[,id], p=bs[,pos], .combine=rbind) %do% {
      dists=cs@biases[name==n & id!=i,.(name,id,distance=abs(pos-p))]
      if (cs@settings$circularize>0)  dists[,distance:=pmin(distance,cs@settings$circularize+1-distance)]
      dists[,dbin:=cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
      dists[!is.na(dbin),.(name=n,id=i,ncross=.N),by=dbin]
    }
    dists
  }
  setkey(zeros,id,name,dbin)
  stopifnot(zeros[,uniqueN(id)]==cs@biases[,.N])
  zeros=rbind(cts[cat=="contact L"][zeros,.(name,id,dbin,cat="contact L",nnz,nposs=2*ncross)],
              cts[cat=="contact R"][zeros,.(name,id,dbin,cat="contact R",nnz,nposs=2*ncross)])
  zeros[is.na(nnz),nnz:=0]
  setkey(zeros,id,name,dbin,cat)
  zeros[,nzero:=nposs-nnz]
  #add positions
  zeros=cs@biases[,.(name,id,pos)][zeros]
  setkeyv(zeros,c("id","name","dbin","cat"))
  #Nkz=csnorm:::get_nzeros_per_decay(cs,ncores)
  #all(zeros[,.(ncross=sum(nposs)/8,nnz=sum(nnz)/2,nzero=sum(nzero)/2),keyby=c("name","dbin")]==Nkz)
  #zbias=csnorm:::get_nzeros_per_cutsite(cs,ncores)
  #all(zeros[,.(nnz=sum(nnz),nzero=sum(nzero)),keyby=c("name","id","cat")]==zbias[,.(name,id,cat,nnz,nzero)])
  return(zeros)
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


#' Prepare for concurrent signal estimation 
#' @keywords internal
prepare_signal = function(cs, base.res) {
  ### build matrix
  #create an empty matrix containing all cells, even those with no cut-site intersection
  signal.bins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+base.res],base.res)
  cs@settings$sbins=signal.bins
  signal.bins=unique(cut(c(signal.bins,head(signal.bins,n=-1)+base.res/2), signal.bins,
                         ordered_result=T, right=F, include.lowest=T,dig.lab=12))
  signal.mat=CJ(name=cs@experiments[,name],bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
  signal.mat[,phi:=0]
  ### build optimization trails
  trails = csnorm:::gfl_compute_trails(signal.mat[,nlevels(bin1)])
  stopifnot(all(signal.mat[,.N,by=name]$N==signal.mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  stopifnot(all(length(V(trails$graph))==signal.mat[,.N,by=name]$N))
  cs@par$signal=signal.mat
  cs@settings$trails=trails
  return(cs)
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
  vars=foreach(var=c("eC","eRJ","eDE","alpha","lambda_iota","lambda_rho","lambda_diag", "lambda1", "lambda2", "eCprime"),
               trans=(c(NA,NA,NA,NA,"log10","log10","log10","log10","log10",NA)),.combine=rbind) %do% get_all_values(cs,var,trans)
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
#' @param init boolean (default is FALSE). Whether to continue an existing
#'  normalization or to start anew (default).
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
run_gauss = function(cs, restart=F, bf_per_kb=30, bf_per_decade=20, bins_per_bf=10, base.res=10000,
                     ngibbs = 20, iter=10000, fit.decay=T, fit.genomic=T, fit.signal=T, fit.disp=T,
                     verbose=T, ncounts=1000000, init_alpha=1e-7, init.dispersion=10,
                     tol.obj=1e-1, type="perf", ncores=1) {
  #basic checks
  type=match.arg(type,c("outer","perf"))
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #
  if (verbose==T) cat("Normalization with fast approximation and",type,"iteration\n")
  setkey(cs@biases, id, name)
  setkey(cs@counts, id1, id2, name)
  if (restart==F) {
    #fresh start
    cs@par=list() #in case we have a weird object
    cs@groups=list()
    #add settings
    cs@settings = c(cs@settings[c("circularize","dmin","dmax","qmin","qmax")],
                    list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf, base.res=base.res,
                         iter=iter, init_alpha=init_alpha, init.dispersion=init.dispersion, tol.obj=tol.obj))
    cs@settings$Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    cs@settings$Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
    #initial guess
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
    #get number of zeros along cut sites and decay
    if(verbose==T) cat("Counting zeros\n")
    stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
    cs@settings$dbins=10**seq(log10(cs@settings$dmin),log10(cs@settings$dmax)+stepsz,stepsz)
    cs@zeros = list(bg=csnorm:::get_nzeros_normalization(cs, ncores=ncores), sig=data.table())
  } else {
    if (verbose==T) cat("Continuing already started normalization with its original settings\n")
    laststep = cs@diagnostics$params[,max(step)]
    init.mean="mean"
  }
  #
  #prepare signal matrix and trails
  if (fit.signal==T & cs@zeros$sig[,.N]==0) { #either restart==F or restart==T but was run with fit.signal==F
    if(verbose==T) cat("Preparing for signal estimation\n")
    cs = csnorm:::prepare_signal(cs, base.res)
    cs@zeros=modifyList(list(sig=csnorm:::get_nzeros_binning(cs, base.res, ncores = ncores)), cs@zeros)
  }
  #
  if(verbose==T) cat("Subsampling counts for dispersion\n")
  subcounts = csnorm:::subsample_counts(cs, ncounts)
  subcounts.weight = merge(cs@zeros$bg[,.(nc=sum(nposs)/8),by=name],subcounts[,.(ns=.N),keyby=name])[,.(name,wt=nc/ns)]
  #gibbs sampling
  for (i in (laststep + 1:ngibbs)) {
    #fit diagonal decay given iota and rho
    if (fit.decay==T) {
      if (verbose==T) cat("Gibbs",i,": Decay ")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_decay(cs, init.mean=init.mean, type=type)))
      if(length(output) == 0) { output = 'ok' }
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="decay", out=output, runtime=a[1]+a[4], type=type)
      if (verbose==T) cat("log-likelihood = ",cs@par$value, "\n")
    }
    #fit iota and rho given diagonal decay
    if (fit.genomic==T) {
      if (verbose==T) cat("Gibbs",i,": Genomic ")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_genomic(cs, init.mean=init.mean, type=type)))
      if(length(output) == 0) { output = 'ok' }
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="bias", out=output, runtime=a[1]+a[4], type=type)
      if (verbose==T) cat("log-likelihood = ",cs@par$value, "\n")
    }
    init.mean="mean"
    #fit signal using sparse fused lasso
    if (fit.signal==T) {
      if (verbose==T) cat("Gibbs",i,": Signal ")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_signal(cs, verbose=verbose, ncores=ncores)))
      if(length(output) == 0) { output = 'ok' }
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="signal", out=output, runtime=a[1]+a[4], type="perf")
      if (verbose==T) cat("BIC = ",cs@par$value, "\n")
    }
    #fit dispersion
    if (fit.disp==T) {
      if (verbose==T) cat("Gibbs",i,": Remaining parameters ")
      a=system.time(output <- capture.output(cs <- csnorm:::csnorm_gauss_dispersion(cs, counts=subcounts, weight=subcounts.weight,
                                                                                    init_alpha=init_alpha, type=type, ncores = ncores)))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="disp", out=output, runtime=a[1]+a[4], type=type)
      if (verbose==T) cat("log-likelihood = ",cs@par$value,"\n")
    }
    if (verbose==T) cat("Gibbs",i,": dispersion = ",cs@par$alpha, " lambda_iota = ",cs@par$lambda_iota, "\n")
    if (i-laststep>1) if (has_converged(cs)) break
  }
  if (verbose==T) cat("Done\n")
  return(cs)
}


