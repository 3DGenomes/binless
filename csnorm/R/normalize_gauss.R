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
  zdecay = cs@zeros[,.(nzero=sum(nzero)/2),keyby=c("name","dbin")] #each contact is counted twice
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
csnorm_gauss_common_muhat_mean = function(cs, zeros, sbins) {
  init=cs@par
  ### positive counts
  #compute means
  cpos=copy(cs@counts)
  bsub=cs@biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  cpos=merge(bsub[,.(id1=id,log_iota,log_rho)],cpos,by="id1",all.x=F,all.y=T)
  cpos=merge(bsub[,.(id2=id,log_iota,log_rho)],cpos,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  cpos=merge(cbind(cs@design[,.(name)],eC=init$eC), cpos, by="name",all.x=F,all.y=T)
  cpos[,c("bin1","bin2","dbin"):=
         list(cut(pos1, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(pos2, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cpos=merge(cpos,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  cpos[,log_mu.base:=eC + log_decay]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  cpos[,log_mu.base:=NULL]
  cpos=rbind(cpos[contact.close>0,.(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd",
                                    count=contact.close, lmu=lmu.close, weight=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.far>0,  .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd",
                                    count=contact.far,   lmu=lmu.far,   weight=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.down>0, .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd",
                                    count=contact.down,  lmu=lmu.down,  weight=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.up>0,   .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd",
                                    count=contact.up,    lmu=lmu.up,    weight=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.far>0,  .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev",
                                    count=contact.far,   lmu=lmu.far,   weight=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.close>0,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev",
                                    count=contact.close, lmu=lmu.close, weight=1, eC, log_decay, log_bias=log_iota2)],
             cpos[contact.down>0, .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev",
                                    count=contact.down,  lmu=lmu.down,  weight=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.up>0,   .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev",
                                    count=contact.up,    lmu=lmu.up,    weight=1, eC, log_decay, log_bias=log_iota2)])
  ### zero counts
  czero = merge(init$decay[,.(name,dbin,log_decay)], zeros, by=c("name","dbin"))
  stopifnot(czero[is.na(log_decay),.N]==0)
  czero = merge(czero[nzero>0],bsub, by.x="id1",by.y="id",all.x=T,all.y=F)
  czero[,log_bias:=ifelse(cat=="contact L", log_iota, log_rho)]
  czero = merge(cbind(cs@design[,.(name)],eC=init$eC), czero, by="name",all.x=F,all.y=T)
  czero[,lmu:=eC + log_decay + log_bias] #we don't average over j
  czero = czero[,.(name,id1,pos1,bin1,bin2,dbin,cat,dir,count=0,lmu,weight=nzero,log_bias,log_decay,eC)]
  cts=rbind(cpos,czero)
  ### add signal
  if (cs@par$signal[,.N]>0) {
    signal=rbind(cs@par$signal,cs@par$signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
    cts=signal[cts,,on=c("name","bin1","bin2")]
    cts[,mu:=exp(lmu+phi)]
  } else {
    cts[,mu:=exp(lmu)]
  }
  ### finalize
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  #ggplot(cts[name=="T47D es 60 MboI 1"&cat=="contact R"])+geom_line(aes(dbin,log_decay,colour=count>0,group=count>0))
  #cts[,c("count","mu","pos","lmu"):=NULL] #not needed downstream
  return(cts)
}

#' Compute decay etahat and weights using previous mean
#' @keywords internal
#' 
csnorm_gauss_decay_muhat_mean = function(cs) {
  cts = csnorm:::csnorm_gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
  cts[,kappaij:=eC+log_decay]
  dbins=cs@settings$dbins
  csd = cts[,.(distance=sqrt(dbins[unclass(dbin)+1]*dbins[unclass(dbin)]),
               kappahat=weighted.mean(z+kappaij, weight/var),
               std=1/sqrt(sum(weight/(2*var))), weight=sum(weight)/2), keyby=c("name", "dbin")] #each count appears twice
  return(csd)
}

#' Generate cubic spline
#' @keywords internal
#' 
csnorm_generate_cubic_spline = function(cutsites, Krow, sparse=F) {
  splinedegree=3 #Cubic spline
  dx = 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree)
  t = min(cutsites) - dx*0.01 + dx * seq(-splinedegree,Krow-splinedegree+3)
  return(spline.des(cutsites, knots = t, outer.ok = T, sparse=sparse)$design)
}
  
#' Optimize decay parameters
#' @keywords internal
#' 
csnorm_gauss_decay_optimize = function(csd, design, Kdiag, original_lambda_diag,
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
    X = csnorm_generate_cubic_spline(log(cutsites), Kdiag, sparse=F)
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
        nBsp = csnorm_generate_cubic_spline(log(ncutsites), Kdiag, sparse=F)
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
      
      nlambda_diag = (Kdiag - 2)/((Kdiag^2)*crossprod(D%*%beta)+1)
      nlambda_diag = sqrt(as.numeric(nlambda_diag))
      
      epsilon = abs(lambda_diag-nlambda_diag)
      lambda_diag = nlambda_diag
      maxiter = maxiter+1
    }
    if (verbose==T) cat("   step",maxiter-1,": lambda_diag",nlambda_diag,"\n")
    
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
    
    decay_out$value = sum(dnorm(kappa_hat, mean = as.array(decay_out$log_mean_counts), sd = as.array(sdl), log = TRUE))
    
  }
  decay_out$eC=as.array(decay_out$eC)
  decay_out$beta_diag=as.matrix(decay_out$beta_diag)
  decay_out$beta_diag_diff=as.matrix(decay_out$beta_diag_diff)
  #make decay data table, reused at next call
  dmat=csd[,.(name,dbin,distance,kappahat,std,ncounts=weight,kappa=decay_out$log_mean_counts,log_decay=decay_out$log_decay)]
  setkey(dmat,name,dbin)
  decay_out$decay=dmat 
  return(decay_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @keywords internal
#' 
csnorm_gauss_decay = function(cs, verbose=T, init.mean="mean", update.eC=T) {
  if (verbose==T) cat(" Decay\n")
  if (init.mean=="mean") {
    csd = csnorm:::csnorm_gauss_decay_muhat_mean(cs)
  } else {
    csd = csnorm:::csnorm_gauss_decay_muhat_data(cs)
  }
  #run optimization
  op = csnorm:::csnorm_gauss_decay_optimize(csd, cs@design, cs@settings$Kdiag, cs@par$lambda_diag,
                                            verbose=verbose, max_perf_iteration=cs@settings$iter,
                                            convergence_epsilon=cs@settings$tol.leg)
  #update par slot
  if (update.eC==F) op$eC=NULL
  cs@par=modifyList(cs@par, op)
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
  zbias = cs@zeros[,.(id=id1,name,cat,nzero)][,.(nzero=sum(nzero)),keyby=c("id","name","cat")][cs@biases[,.(name,id,pos)]]
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
  cts.common = csnorm:::csnorm_gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
  cts.common[,etaij:=eC+log_bias]
  cts = cts.common[,.(etahat=weighted.mean(z+etaij, weight/var),std=sqrt(2/sum(weight/var))),
                   keyby=c("id1","name","cat")]
  cts = merge(cs@biases[,.(name,id,pos)],cts,by.x=c("name","id"),by.y=c("name","id1"))
  setkey(cts, id, name, cat)
  return(list(bts=bts,cts=cts))
}


csnorm_gauss_genomic_optimize = function(bts, cts, biases, design, Krow, sbins,
                                         lambda_iota, lambda_rho, verbose=T,
                                         max_perf_iteration=1000, convergence_epsilon=1e-5) {
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
    Bsp = csnorm_generate_cubic_spline(cutsites, Krow, sparse=T)
    X = rbind(cbind(Bsp/2,Bsp/2),bdiag(Bsp,Bsp),bdiag(Bsp,Bsp))
    if (length(sbins)<=2) {
      centering=Matrix(rep(1,SD),ncol=1)
    } else {
      sbinned=biases[bbegin[1]:(bbegin[2]-1), cut(pos, sbins, ordered_result=T,
                                                  right=F, include.lowest=T,dig.lab=12)]
      centering=Matrix(model.matrix(~ 0+sbinned))
      stopifnot(dim(centering)==c(SD,length(sbins)-1))
    }
    W=cbind(rbind(Matrix(0,nrow=3*SD,ncol=ncol(centering)),centering,Matrix(0,nrow=SD,ncol=ncol(centering))),
            rbind(Matrix(0,nrow=4*SD,ncol=ncol(centering)),centering))
    diags = list(rep(1,Krow), rep(-2,Krow))
    D1 = bandSparse(Krow-2, Krow, k=0:2, diagonals=c(diags, diags[1]))
    if (!is.null(lambda_iota) && !is.null(lambda_rho)) {
      lambda_iota = lambda_iota[uXB]
      lambda_rho = lambda_rho[uXB]
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
        nBsp  = csnorm_generate_cubic_spline(ncutsites, Krow, sparse=T)
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
    
    epsilon = 1
    maxiter = 0
    while(epsilon > convergence_epsilon && maxiter < max_perf_iteration) {
      D = bdiag(lambda_iota*D1,lambda_rho*D1)
      DtD = crossprod(D)
      if (maxiter==0) {
        cholA = Cholesky(tmp_X_S_m2_X + Krow^2*DtD, LDL=F, super=NA) 
        stopifnot(!isLDL(cholA)) #do LLt cholesky so we can use crossprod for Gamma_v
      } else {
        cholA = update(cholA,tmp_X_S_m2_X + Krow^2*DtD)
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
      
      beta = beta_y - beta_U%*%e
      beta_iota = beta[1:Krow]
      beta_rho = beta[(Krow+1):(2*Krow)]
      
      nlambda_iota = (Krow - 2)/((Krow**2)*crossprod(D1%*%beta_iota)+1e6)
      nlambda_iota = sqrt(as.numeric(nlambda_iota))
      nlambda_rho = (Krow - 2)/((Krow**2)*crossprod(D1%*%beta_rho)+1e6)
      nlambda_rho = sqrt(as.numeric(nlambda_rho))
      
      epsilon = max(abs(lambda_iota-nlambda_iota),abs(lambda_rho-nlambda_rho))
      lambda_iota = nlambda_iota
      lambda_rho = nlambda_rho
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
    genomic_out$lambda_iota = c(genomic_out$lambda_iota,lambda_iota)
    genomic_out$lambda_rho = c(genomic_out$lambda_rho,lambda_rho)
    genomic_out$beta_iota_diff = as.matrix(genomic_out$beta_iota_diff)
    genomic_out$beta_rho_diff = as.matrix(genomic_out$beta_rho_diff)
    
    means = cbind(genomic_out$log_mean_RJ,genomic_out$log_mean_DL,genomic_out$log_mean_DR,genomic_out$log_mean_cleft,genomic_out$log_mean_cright)
    mus = cbind(sd_RJ,sd_DL,sd_DR,sd_L,sd_R)
    genomic_out$value = sum(dnorm(etas, mean = as.array(means), sd = as.array(mus), log = TRUE))
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
  return(genomic_out)
}

#' Single-cpu simplified fitting for iota and rho
#' @param if a single value, use data for estimate of mu and that value as a
#'   dispersion, otherwise it's a list with parameters to compute the mean from
#' @keywords internal
#'   
csnorm_gauss_genomic = function(cs, verbose=T, init.mean="mean", update.exposures=T) {
  if (verbose==T) cat(" Genomic\n")
  if (init.mean=="mean") {
    a = csnorm:::csnorm_gauss_genomic_muhat_mean(cs)
  } else {
    a = csnorm:::csnorm_gauss_genomic_muhat_data(cs)
  }
  #run optimization
  op = csnorm:::csnorm_gauss_genomic_optimize(a$bts, a$cts, cs@biases, cs@design, cs@settings$Krow, cs@settings$sbins,
                                              cs@par$lambda_iota, cs@par$lambda_rho, verbose=verbose,
                                              max_perf_iteration=cs@settings$iter,
                                              convergence_epsilon=cs@settings$tol.leg)
  #update par slot
  if (update.exposures==F) {
    op$eC=NULL
    op$eRJ=NULL
    op$eDE=NULL
  }
  cs@par=modifyList(cs@par, op)
  return(cs)
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
csnorm_gauss_dispersion = function(cs, counts, weight=cs@design[,.(name,wt=1)], verbose=T, ncores=ncores) {
  if (verbose==T) cat(" Dispersion\n")
  #predict all means and put into table
  counts=csnorm:::csnorm_predict_all_parallel(cs,counts,ncores = ncores)
  setkeyv(counts,key(cs@counts))
  stopifnot(cs@biases[,.N]==length(cs@par$log_iota))
  #add signal contribution if available
  if (length(cs@settings$sbins)>2) {
    counts[,bin1:=cut(pos1, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    counts[,bin2:=cut(pos2, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    counts=cs@par$signal[counts,,on=key(cs@par$signal)]
    counts[,c("log_mean_cclose","log_mean_cfar","log_mean_cup","log_mean_cdown"):=
             list(log_mean_cclose+phi,log_mean_cfar+phi,log_mean_cup+phi,log_mean_cdown+phi)]
  }
  if (verbose==T) cat("  predict\n")
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
  out=capture.output(op<-optimize_stan_model(model=csnorm:::stanmodels$gauss_dispersion,
                                             data=data, iter=cs@settings$iter, tol_obj=cs@settings$tol.leg,
                                             verbose=verbose, init=init, init_alpha=cs@settings$init_alpha))
  cs@par=modifyList(cs@par, op$par[c("eRJ","eDE","alpha")])
  #cs@par$eC=cs@par$eC+op$par$eC_sup
  if (verbose==T) cat("  fit: dispersion",cs@par$alpha,"\n")
  #
  #compute log-posterior
  Krow=cs@settings$Krow
  Kdiag=cs@settings$Kdiag
  cs@par$value = op$value + (Krow-2)/2*sum(log(cs@par$lambda_iota/exp(1))+log(cs@par$lambda_rho/exp(1))) +
    (Kdiag-2)/2*sum(log(cs@par$lambda_diag/exp(1)))
  return(cs)
}

#' fit signal using sparse fused lasso
#' @keywords internal
#' 
csnorm_gauss_signal = function(cs, verbose=T, ncores=ncores) {
  if (verbose==T) cat(" Signal\n")
  cts = csnorm:::csnorm_gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
  mat = cts[,.(phihat=weighted.mean(z+phi, weight/var),
               phihat.var=2/sum(weight/var),
               ncounts=sum(weight)),keyby=c("name","bin1","bin2")][cs@par$signal[,.(name,bin1,bin2)],,on=c("name","bin1","bin2")]
  mat[is.na(phihat),c("phihat","phihat.var","ncounts"):=list(1,Inf,0)] #bins with no detectable counts
  mat[,c("valuehat","weight"):=list(phihat,1/phihat.var)]
  setkey(mat,name,bin1,bin2)
  #
  if (cs@par$signal[bin1==bin2,any(phi!=0)]) { #do not remove anything at first iteration
    if (verbose==T) cat("  Remove bad bins\n")
    tmp=cs@par$signal[bin1==bin2,.(is.bad=all(phi==0)),by=bin1][is.bad==T,unclass(bin1)]
    if (length(cs@par$badbins)>0) tmp=union(tmp,cs@par$badbins)
    if (length(tmp)>cs@settings$qmax*cs@par$signal[bin1==bin2&name==name[1],.N]) {
      cat("Refusing to remove more signal rows than the fraction qmax=",cs@settings$qmax,"allows. Attempted:",
          length(tmp),"allowed:",as.integer(cs@settings$qmax*cs@par$signal[bin1==bin2&name==name[1],.N]),"\n")
    } else {
      cs@par$badbins=tmp
    }
    if(length(cs@par$badbins)>0) {
      cat("removing",length(cs@par$badbins),"bad signal rows\n")
      mat[unclass(bin1) %in% cs@par$badbins | unclass(bin2) %in% cs@par$badbins, weight:=0]
    }
  }
  #
  if (verbose==T) cat("  predict\n")
  #ggplot(mat[bin1==bin2])+geom_point(aes(bin1,phi,colour=bin1%in%cs@par$badbins))+facet_wrap(~name)
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phihat))+facet_wrap(~name)
  #ggplot(mat[unclass(bin2)<23&unclass(bin1)>13])+geom_raster(aes(bin1,bin2,fill=phihat))+#geom_raster(aes(bin2,bin1,fill=phi))+
  #  facet_wrap(~name)+scale_fill_gradient2()
  #ggplot(mat[unclass(bin2)<23&unclass(bin1)>13])+geom_raster(aes(bin1,bin2,fill=weight))+
  #  facet_wrap(~name)+scale_fill_gradient2()
  #mat[unclass(bin2)<23&unclass(bin1)>13][bin1==bin2]
  #ggplot(cts.cp)+geom_violin(aes(bin1,count))
  #
  #perform fused lasso on signal
  groupnames=mat[,unique(name)]
  registerDoParallel(cores=ncores)
  params = foreach(g=groupnames, .combine=rbind) %dopar%
    csnorm:::csnorm_fused_lasso(mat[name==g], cs@settings$trails, fixed=F, positive=T,
                                tol=cs@settings$tol.leg, ncores=ncores, verbose=verbose)
  #compute matrix at new params
  #save(mat,params,file=paste0("mat_step_",step,".RData"))
  mat = foreach (g=groupnames, .combine=rbind) %dopar% {
    p=params[name==g]
    matg=mat[name==g]
    matg[,value:=csnorm:::gfl_get_value(valuehat, weight, cs@settings$trails, p$lambda1, p$lambda2, p$eCprime)]
    matg
  }
  #store new signal in cs and update eC
  #mat[,phi.old:=phi]
  mat[,phi:=value]
  #ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+
  #  scale_fill_gradient2()
  #ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=phi==0))
  setkey(mat,name,bin1,bin2)
  cs@par$signal=mat[,.(name,bin1,bin2,phi)]
  params=merge(cbind(cs@design[,.(name)],eC=cs@par$eC), params, by="name",all=T)
  cs@par$eC=as.array(params[,eC+eCprime])
  cs@par$eCprime=as.array(params[,eCprime])
  cs@par$lambda1=as.array(params[,lambda1])
  cs@par$lambda2=as.array(params[,lambda2])
  cs@par$value = params[,sum(BIC)]
  if (verbose==T) cat("  fit: lambda1",cs@par$lambda1[1],"lambda2",cs@par$lambda2[1],"\n")
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

#' count number of zeros in a given cut site, distance bin and signal bin
#' @keywords internal
#' 
get_nzeros = function(cs, sbins, ncores=1) {
  stopifnot(cs@counts[id1>=id2,.N]==0)
  #positive counts: group per cut site and signal / distance bin
  cts=melt(cs@counts[,.(name,pos1,pos2,distance,contact.close,contact.down,contact.far,contact.up)],
           id.vars=c("name","pos1","pos2","distance"))[value>0]
  cts[,c("bin1","bin2","dbin"):=
        list(cut(pos1, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
             cut(pos2, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
             cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cts=rbind(cts[,.(name,pos1,bin1,bin2,dbin,variable)][
              ,.(nnz=.N,dir="fwd"),keyby=c("name","pos1","bin1","bin2","dbin","variable")],
            cts[,.(name,pos1=pos2,bin1=bin2,bin2=bin1,dbin,variable)][
              ,.(nnz=.N,dir="rev"),keyby=c("name","pos1","bin1","bin2","dbin","variable")])
  cts[,cat:=ifelse(dir=="fwd",ifelse(variable %in% c("contact.up","contact.far"), "contact L", "contact R"),
                              ifelse(variable %in% c("contact.up","contact.close"), "contact L", "contact R"))]
  cts=cts[,.(nnz=sum(nnz)),keyby=c("name","pos1","bin1","bin2","dbin","cat","dir")]
  #Count the number of crossings per distance bin
  #looping over IDs avoids building NxN matrix
  biases=cs@biases[,.(name,id,pos)]
  biases[,bin:=cut(pos, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
  registerDoParallel(cores=ncores)
  chunksize=cs@biases[,ceiling(.N/(10*ncores))]
  nchunks=cs@biases[,ceiling(.N/chunksize)]
  crossings = foreach(chunk=1:nchunks, .combine=rbind) %dopar% {
    bs=biases[((chunk-1)*chunksize+1):min(.N,chunk*chunksize)]
    foreach(i=bs[,id], n=bs[,name], p=bs[,pos], b=bs[,bin], .combine=rbind) %do% {
      crossings = biases[name==n&pos!=p,.(name,pos2=pos,bin2=bin,distance=abs(pos-p))]
      if (cs@settings$circularize>0)  crossings[,distance:=pmin(distance,cs@settings$circularize+1-distance)]
      crossings[,c("dbin","dir"):=list(
        cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12),
        ifelse(pos2>p,"fwd","rev"))]
      crossings[distance>=cs@settings$dmin,.(id1=i, pos1=p,bin1=b,ncross=.N),by=c("name","bin2","dbin","dir")]
    }
  }
  #now merge with positive counts and deduce number of zeros
  zeros = rbind(
    merge(crossings,cts[cat=="contact L"],by=c("name","pos1","bin1","bin2","dbin","dir"),all=T)[
      ,.(name,id1,pos1,bin1,bin2,dbin,dir,cat="contact L",ncross,nnz)],
    merge(crossings,cts[cat=="contact R"],by=c("name","pos1","bin1","bin2","dbin","dir"),all=T)[
      ,.(name,id1,pos1,bin1,bin2,dbin,dir,cat="contact R",ncross,nnz)])
  zeros[is.na(nnz),nnz:=0]
  zeros[,nzero:=2*ncross-nnz]
  stopifnot(zeros[is.na(ncross),.N==0])
  stopifnot(zeros[nzero<0,.N==0])
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
prepare_signal_estimation = function(biases, names, base.res) {
  ### build matrix
  #create an empty matrix containing all cells, even those with no cut-site intersection
  sbins=seq(biases[,min(pos)-1],biases[,max(pos)+1+base.res],base.res)
  signal.bins=unique(cut(c(sbins,head(sbins,n=-1)+base.res/2), sbins,
                         ordered_result=T, right=F, include.lowest=T,dig.lab=12))
  signal.mat=CJ(name=names,bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
  signal.mat[,phi:=0]
  setkey(signal.mat,name,bin1,bin2)
  ### build optimization trails
  trails = csnorm:::gfl_compute_trails(signal.mat[,nlevels(bin1)])
  stopifnot(all(signal.mat[,.N,by=name]$N==signal.mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  stopifnot(all(length(V(trails$graph))==signal.mat[,.N,by=name]$N))
  return(list(signal=signal.mat,trails=trails,sbins=sbins))
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
plot_diagnostics = function(cs, start=1) {
  plot=ggplot(cs@diagnostics$params[step>=start,.(step,leg,value)])+
    geom_line(aes(step,value))+geom_point(aes(step,value))+facet_wrap(~leg, scales = "free")+
    theme(legend.position="bottom")
  vars=foreach(var=c("eC","eRJ","eDE","alpha","lambda_iota","lambda_rho","lambda_diag", "lambda1", "lambda2", "eCprime"),
               trans=(c(NA,NA,NA,NA,"log10","log10","log10","log10","log10",NA)),.combine=rbind) %do% get_all_values(cs,var,trans)
  plot2=ggplot(vars[step>=start])+geom_line(aes(step,value))+
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
#' @param iter positive integer. Maximum number of optimization steps for each leg.
#' @param fit.decay,fit.genomic,fit.signal,fit.disp boolean. Set to FALSE only for diagnostics.
#' @param verbose Display progress if TRUE
#' @param init_alpha positive numeric, default 1e-5. Initial step size of LBFGS
#'   line search (dispersion leg).
#' @param init.dispersion positive numeric. Value of the dispersion to use initially.
#' @param tol.obj positive numeric (default 1e-1). Convergence tolerance on changes in the four likelihoods.
#' @param tol.val positive numeric (default 1e-3). Convergence tolerance on (relative) changes in the parameters.
#' @param tol.leg positive numeric (default 1e-5). Convergence tolerance on changes in the likelihood within each leg.
#' @param ncores positive integer (default 1). Number of cores to use.
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
#' 
run_gauss = function(cs, restart=F, bf_per_kb=30, bf_per_decade=20, bins_per_bf=10, base.res=10000,
                     ngibbs = 20, iter=1000, fit.decay=T, fit.genomic=T, fit.signal=T, fit.disp=T,
                     verbose=T, ncounts=1000000, init_alpha=1e-7, init.dispersion=10,
                     tol.obj=1e-1, tol.val=1e-3, tol.leg=1e-5, ncores=1) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #
  if (verbose==T) cat("Normalization with fast approximation and performance iteration\n")
  setkey(cs@biases, id, name)
  setkey(cs@counts, id1, id2, name)
  if (restart==F) {
    #fresh start
    cs@par=list() #in case we have a weird object
    cs@groups=list()
    #add settings
    cs@settings = c(cs@settings[c("circularize","dmin","dmax","qmin","qmax")],
                    list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf, base.res=base.res,
                         iter=iter, init_alpha=init_alpha, init.dispersion=init.dispersion, tol.obj=tol.obj,
                         tol.val=tol.val, tol.leg=tol.leg))
    cs@settings$Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    cs@settings$Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
    stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
    cs@settings$dbins=10**seq(log10(cs@settings$dmin),log10(cs@settings$dmax)+stepsz,stepsz)
    #initial guess
    if (verbose==T) cat("No initial guess provided\n")
    cs@diagnostics=list()
    laststep=0
    init.mean="data"
    cs=csnorm:::fill_parameters_perf(cs, dispersion=init.dispersion, fit.decay=fit.decay,
                                     fit.genomic=fit.genomic, fit.disp=fit.disp)
    #prepare signal matrix and trails
    if (fit.signal==T) {
      if(verbose==T) cat("Preparing for signal estimation\n")
      stuff = csnorm:::prepare_signal_estimation(cs@biases, cs@experiments[,name], base.res)
      cs@par$signal=stuff$signal
      cs@settings$sbins=stuff$sbins
      cs@settings$trails=stuff$trails
    } else {
      cs@settings$sbins=cs@biases[,c(min(pos)-1,max(pos)+1)]
      cs@par$signal=cs@biases[,.(phi=0,bin1=cut(pos[1], cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                                 bin2=cut(pos[1], cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)),by=name]
      setkey(cs@par$signal,name,bin1,bin2)
      cs@par$settings$trails=NA
    }
    #get number of zeros along cut sites and decay
    if(verbose==T) cat("Counting zeros\n")
    cs@zeros = csnorm:::get_nzeros(cs, cs@settings$sbins, ncores=ncores)
    #update eC everywhere in the first step
    update.exposures=T
  } else {
    if (verbose==T) cat("Continuing already started normalization with its original settings\n")
    laststep = cs@diagnostics$params[,max(step)]
    init.mean="mean"
    update.exposures=(fit.signal==T)
  }
  #
  if(verbose==T) cat("Subsampling counts for dispersion\n")
  subcounts = csnorm:::subsample_counts(cs, ncounts)
  subcounts.weight = merge(cs@zeros[,.(nc=sum(ncross)/4),by=name],subcounts[,.(ns=.N),keyby=name],by="name")[,.(name,wt=nc/ns)]
  #gibbs sampling
  for (i in (laststep + 1:ngibbs)) {
    if (verbose==T) cat("\n### Iteration",i,"\n")
    #fit diagonal decay given iota and rho
    if (fit.decay==T) {
      a=system.time(cs <- csnorm:::csnorm_gauss_decay(cs, init.mean=init.mean, update.eC=update.exposures))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="decay", runtime=a[1]+a[4])
      if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
    }
    #fit iota and rho given diagonal decay
    if (fit.genomic==T) {
      a=system.time(cs <- csnorm:::csnorm_gauss_genomic(cs, init.mean=init.mean, update.exposures=update.exposures))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="bias", runtime=a[1]+a[4])
      if (verbose==T) cat("  log-likelihood = ",cs@par$value, "\n")
    }
    init.mean="mean"
    #fit dispersion
    if (fit.disp==T) {
      a=system.time(cs <- csnorm:::csnorm_gauss_dispersion(cs, counts=subcounts, weight=subcounts.weight,
                                                           ncores = ncores))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="disp", runtime=a[1]+a[4])
      if (verbose==T) cat("  log-likelihood = ",cs@par$value,"\n")
    }
    #fit signal using sparse fused lasso
    if (fit.signal==T) {
      update.exposures=F
      a=system.time(cs <- csnorm:::csnorm_gauss_signal(cs, verbose=verbose, ncores=ncores))
      cs@diagnostics$params = csnorm:::update_diagnostics(cs, step=i, leg="signal", runtime=a[1]+a[4])
      if (verbose==T) cat("  BIC = ",cs@par$value, "\n")
    }
    #check for convergence
    if (i>1) if (has_converged(cs)) {
      if (verbose==T) cat("Normalization has converged\n")
      break
    }
  }
  if (verbose==T) cat("Done\n")
  return(cs)
}


