#' @include binless.R
NULL

#' Compute means for positive and zero counts using previous params
#' 
#' @return data.table containing the following columns
#' - name: name of the dataset
#' - id1, pos1, bin1: coordinates of the cut site
#' - bin2, dbin: signal/distance bin in which we looked at the intersections
#' - dir: fwd (rev) for contacts with downstream (resp upstream) cut-sites
#' - cat: whether we consider contacts on the left (contact L) or on the right (contact R) of this cut site
#' - count: the value of the count in this signal/distance/direction/category bin.
#'  We discard anything below cs@settings$dmin. Note that sum(count) is twice the total number of counts (per dataset).
#' - lmu.nosig: the log(mean) of the background model (exposure+bias+decay)
#' - phi: the log(signal)
#' - mu: the mean, including signal
#' - weight: how many observations in this signal/distance/direction/category bin have this count.
#'  Note that sum(weight) is twice the total number of observable counts (per dataset).
#' - eC, log_decay, log_bias: exposure, log(decay) and log(bias) of this dataset. The bias is either iota or rho depending on cat.
#' - z: count/mu-1
#' - var: 1/mu+1/dispersion
#'   
#' @keywords internal
#' 
gauss_common_muhat_mean = function(cs, zeros, sbins) {
  init=cs@par
  ### positive counts (twice, in both directions)
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
                                    count=contact.close, lmu.nosig=lmu.close, weight=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.far>0,  .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd",
                                    count=contact.far,   lmu.nosig=lmu.far,   weight=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.down>0, .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd",
                                    count=contact.down,  lmu.nosig=lmu.down,  weight=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.up>0,   .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd",
                                    count=contact.up,    lmu.nosig=lmu.up,    weight=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.far>0,  .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev",
                                    count=contact.far,   lmu.nosig=lmu.far,   weight=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.close>0,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev",
                                    count=contact.close, lmu.nosig=lmu.close, weight=1, eC, log_decay, log_bias=log_iota2)],
             cpos[contact.down>0, .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev",
                                    count=contact.down,  lmu.nosig=lmu.down,  weight=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.up>0,   .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev",
                                    count=contact.up,    lmu.nosig=lmu.up,    weight=1, eC, log_decay, log_bias=log_iota2)])
  ### zero counts (twice, in both directions)
  czero = init$decay[zeros,.(name,dbin,bin1,bin2,cat,dir,id1,pos1,weight=nzero,log_decay),on=c("name","dbin")]
  setkeyv(czero,key(zeros))
  stopifnot(czero[is.na(log_decay),.N]==0)
  setnames(bsub,"id","id1")
  czero = bsub[czero,on="id1"]
  czero[,log_bias:=ifelse(cat=="contact L", log_iota, log_rho)]
  czero[,c("log_iota","log_rho"):=NULL]
  czero = cbind(cs@design[,.(name)],eC=init$eC)[czero,on="name"]
  czero[,lmu.nosig:=eC + log_decay + log_bias] #we don't average over j
  czero[,count:=0]
  cts=rbind(cpos,czero)
  setkeyv(cts,key(zeros))
  ### add signal
  signal = binless:::get_signal_matrix(cs, resolution = sbins[2]-sbins[1], groups=cs@experiments[,.(name,groupname=name)])
  signal=rbind(signal[,.(name,bin1,bin2,phi)],signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
  cts=signal[cts,,on=c("name","bin1","bin2")]
  cts[,mu:=exp(lmu.nosig+phi)]
  ### finalize
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  setkey(cts,name,id1,pos1,bin1,bin2,dbin,dir,cat)
  #ggplot(cts[name=="T47D es 60 MboI 1"&cat=="contact R"])+geom_line(aes(dbin,log_decay,colour=count>0,group=count>0))
  return(cts)
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

#' Generate cubic spline
#' @keywords internal
#' 
generate_cubic_spline = function(cutsites, Krow, sparse=F) {
  splinedegree=3 #Cubic spline
  dx = 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree)
  t = min(cutsites) - dx*0.01 + dx * seq(-splinedegree,Krow-splinedegree+3)
  return(splineDesign(cutsites, knots = t, ord = splinedegree + 1, outer.ok = F, sparse=sparse))
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
    X = generate_cubic_spline(log(cutsites), Kdiag, sparse=F)
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
        nBsp = generate_cubic_spline(log(ncutsites), Kdiag, sparse=F)
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
      fit = solve.QP(At.PD, tmp_X_S_m2_k, -Ct, meq = Dsets)
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

#' Compute means for a given counts matrix
#' @keywords internal
compute_means = function(cs, counts) {
  #compute background
  init=cs@par
  cpos=copy(counts)
  bsub=cs@biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  cpos=merge(bsub[,.(id1=id,log_iota,log_rho)],cpos,by="id1",all.x=F,all.y=T)
  cpos=merge(bsub[,.(id2=id,log_iota,log_rho)],cpos,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  cpos=merge(cbind(cs@design[,.(name)],eC=init$eC), cpos, by="name",all.x=F,all.y=T)
  cpos[,c("bin1","bin2","dbin"):=
         list(cut(pos1, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(pos2, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cpos=merge(cpos,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  #compute signal
  if (cs@par$signal[,.N]>0 && length(cs@settings$sbins)>2) {
    signal = binless:::get_signal_matrix(cs, resolution = cs@settings$base.res, groups=cs@experiments[,.(name,groupname=name)])
    signal = rbind(signal[,.(name,bin1,bin2,phi)],signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
    cpos = signal[cpos,,on=c("name","bin1","bin2")]
  } else {
    cpos[,phi:=0]
  }
  #assemble
  cpos[,log_mu.base:=eC + log_decay + phi]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  
  cpos=cpos[,.(id1,id2,name,pos1,pos2,distance,contact.close,contact.far,contact.up,contact.down,
               log_decay,lmu.close,lmu.far,lmu.up,lmu.down)]
  setkeyv(cpos,key(cs@counts))
  cpos
}

#' Single-cpu simplified fitting for exposures and dispersion
#' @keywords internal
#' 
gauss_dispersion = function(cs, counts, weight=cs@design[,.(name,wt=1)], verbose=T) {
  if (verbose==T) cat(" Dispersion\n")
  #predict all means and put into table
  counts = binless:::compute_means(cs,counts)
  stopifnot(cs@biases[,.N]==length(cs@par$log_iota))
  #
  #fit dispersion and exposures
  if (verbose==T) cat("  predict\n")
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
               log_mean_cclose=counts[,lmu.close], log_mean_cfar=counts[,lmu.far],
               log_mean_cup=counts[,lmu.up], log_mean_cdown=counts[,lmu.down])
  init=list(eC_sup=as.array(counts[,log(mean(contact.close/exp(lmu.close))),by=name][,V1]),
            eRJ=as.array(cs@biases[,.(name,frac=rejoined/exp((cs@par$log_iota+cs@par$log_rho)/2))][,log(mean(frac)),by=name][,V1]),
            eDE=as.array(cs@biases[,.(name,frac=(dangling.L/exp(cs@par$log_iota)+dangling.R/exp(cs@par$log_rho))/2)][
              ,log(mean(frac)),by=name][,V1]))
  init$mu=mean(exp(init$eC_sup[1]+counts[name==name[1],lmu.close]))
  init$alpha=max(0.001,1/(var(counts[name==name[1],contact.close]/init$mu)-1/init$mu))
  init$mu=NULL
  out=capture.output(op<-optimize_stan_model(model=binless:::stanmodels$gauss_dispersion, tol_param=cs@par$tol_disp,
                                             data=data, iter=cs@settings$iter, verbose=verbose, init=init,
                                             init_alpha=1e-9))
  #restrict tolerance if needed
  precision = max(abs(c(op$par[c("eRJ","eDE","alpha")],recursive=T) - c(cs@par[c("eRJ","eDE","alpha")],recursive=T)))
  cs@par$tol_disp = min(cs@par$tol_disp, max(cs@settings$tol, precision/10))
  #update parameters
  cs@par=modifyList(cs@par, op$par[c("eRJ","eDE","alpha")])
  #cs@par$eC=cs@par$eC+op$par$eC_sup
  #
  #compute log-posterior
  Krow=cs@settings$Krow
  Kdiag=cs@settings$Kdiag
  cs@par$value = op$value + (Krow-2)/2*sum(log(cs@par$lambda_iota/exp(1))+log(cs@par$lambda_rho/exp(1))) +
    (Kdiag-2)/2*sum(log(cs@par$lambda_diag/exp(1)))
  if (verbose==T) {
    cat("  fit: dispersion",cs@par$alpha,"\n")
    cat("  log-likelihood = ",cs@par$value,"\n")
  }
  return(cs)
}

#' fit signal using sparse fused lasso
#' @keywords internal
#' 
gauss_signal_muhat_mean = function(cs, cts.common) {
  cts = cts.common[,.(name,bin1=pmin(bin1,bin2),bin2=pmax(bin1,bin2),count,z,var,mu,lmu.nosig,
                      log_decay,log_bias,weight=weight/2)]
  stopifnot(cts[,all(bin1<=bin2)])
  return(cts)
}

#get metadata for signal calculation
#consists of outliers, which should be discarded for signal detection,
#and information relative to constraining the signal wrt the decay
#' @keywords internal
#' 
get_signal_metadata = function(cs, cts, resolution) {
  cts_compressed = cts[,.(weight=sum(weight), z=weighted.mean(z,weight/var),var=1/weighted.mean(1/var,weight)),
                       keyby=c("name","bin1","bin2")]
  #biases
  # bad.biases=rbind(cts_compressed[,.(bin1,bin2,weight,var,z)],cts_compressed[,.(bin1=bin2,bin2=bin1,weight,var,z)])[
  #   ,.(z=sum(weight*z/var)/sum(weight^2/var^2)),by=bin1]
  # bad.biases[,z:=scale(z)]
  # bad.biases[,is.out:=-abs(z)<qnorm(cs@settings$qmin)]
  # #ggplot(bad.biases)+geom_point(aes(bin1,z,colour=is.out))
  # bad.rows=bad.biases[is.out==T,bin1]
  # if (bad.biases[,sum(is.out)/.N>0.1]) cat(" Warning: removing ",bad.biases[,100*sum(is.out)/.N],"% of all rows!\n")
  # #decay
  cts_compressed[,diag.idx:=unclass(bin2)-unclass(bin1)]
  # bad.decays=cts_compressed[,.(z=sum(weight*z/var)/sum(weight^2/var^2)),by=diag.idx]
  # bad.decays[,z:=scale(z)]
  # bad.decays[,is.out:=-abs(z)<qnorm(cs@settings$qmin)]
  diag.rm = ceiling(cs@settings$dmin/resolution)
  # bad.decays[diag.idx<=diag.rm,is.out:=T]
  # #ggplot(bad.decays)+geom_point(aes(diag.idx,z,colour=is.out))
  # bad.diagonals=bad.decays[is.out==T,diag.idx]
  # if (bad.decays[,sum(is.out)/.N>0.1]) cat(" Warning: removing ",bad.decays[,100*sum(is.out)/.N],"% of all counter diagonals!\n")
  #orthogonality
  dx = 1.01*(log10(cs@settings$dmax)-log10(cs@settings$dmin))/(cs@settings$Kdiag-3)
  orth=data.table(diag.idx=0:cts_compressed[,max(diag.idx)],
                  segment=ceiling((log10(0:cts_compressed[,max(diag.idx)]*resolution)-log10(cs@settings$dmin))/(dx/2))) #double constraint
  orth[,rank:=frank(segment,ties.method = "dense")-1]
  return(list(bad.diagonals=0:diag.rm, bad.rows=as.integer(), diag.grp=orth[,rank]))
  #return(list(bad.diagonals=bad.diagonals,bad.rows=bad.rows, diag.grp=orth[,rank]))
}


#' make a summarized cts matrix used for one IRLS iteration
#' @keywords internal
#' 
compress_cts = function(cs, cts) {
  cts.new = cts[,.(phihat=weighted.mean(count/mu-1,weight/(2*var)), weight=sum(weight/(2*var)),
                   log_decay=weighted.mean(log_decay,weight)), keyby=c("name","bin1","bin2")]
  cts.new[,count:=phihat+2] #or any strictly positive number (but beware of exp(lmu.nosig) )
  cts.new = cs@par$signal[,.(name,bin1,bin2,phi)][cts.new]
  cts.new[,.(name,bin1,bin2,count,lmu.nosig=log(count/(phihat+1))-phi, log_decay, weight=2*weight*(1/cs@par$alpha + (1+phihat)/count))]
}

#' fit signal using sparse fused lasso
#' @keywords internal
#' 
gauss_signal = function(cs, cts.common, verbose=T, ncores=1, fix.lambda1=F, fix.lambda1.at=NA,
                               fix.lambda2=F, fix.lambda2.at=NA) {
  if (verbose==T) cat(" Signal\n")
  #cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
  cts = binless:::gauss_signal_muhat_mean(cs, cts.common)
  metadata = binless:::get_signal_metadata(cs, cts, cs@settings$base.res)
  cts = binless:::compress_cts(cs, cts)
  #
  if (verbose==T) cat("  predict\n")
  #perform fused lasso on signal
  groupnames=cts[,unique(name)]
  nbins=length(cs@settings$sbins)-1
  csigs = foreach(g=groupnames) %do% {
    csig=new("CSbsig", mat=cs@par$signal[name==g], cts=cts[name==g],
             settings=list(metadata=metadata,
                           nbins=nbins, dispersion=cs@par$alpha,
                           last.beta=cs@par$signal[name==g,beta],
                           tol.val=cs@par$tol_signal, nperf=1)) #only one IRLS iteration
    csig@state = binless:::gfl_compute_initial_state(csig, diff=F)
    csig
  }
  registerDoParallel(cores=min(ncores,length(groupnames)))
  params = foreach(csig=csigs, .combine=rbind, .export=c("verbose","fix.lambda1","fix.lambda1.at",
                                                         "fix.lambda2","fix.lambda2.at")) %dopar% {
    binless:::fused_lasso(csig, positive=T, fixed=F, constrained=T, verbose=verbose,
                                fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                                fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at)
  }
  stopImplicitCluster()
  #compute matrix at new params
  mat = rbindlist(params[,mat])
  #store new signal in cs and update eC
  #ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+
  #  scale_fill_gradient2()+coord_fixed()
  #ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=phi==0))
  setkey(mat,name,bin1,bin2)
  #restrict tolerance if needed
  precision = max(abs(mat[,beta]-cs@par$beta.phi))
  cs@par$tol_signal = min(cs@par$tol_signal, max(cs@settings$tol, precision/10))
  #set new parameters
  cs@par$signal=mat[,.(name,bin1,bin2,phihat,weight,ncounts,phi,phi.unconstr,beta,diag.grp,diag.idx)]
  cs@par$beta.phi=mat[,beta]
  params=merge(cbind(cs@design[,.(name)],eC=cs@par$eC), params, by="name",all=T)
  cs@par$eC=as.array(params[,eC+eCprime])
  cs@par$eCprime=as.array(params[,eCprime])
  cs@par$lambda1=as.array(params[,lambda1])
  cs@par$eCprime.unconstr=as.array(params[,eCprime.unconstr])
  cs@par$lambda1.unconstr=as.array(params[,lambda1.unconstr])
  cs@par$lambda2=as.array(params[,lambda2])
  cs@par$value = params[,sum(BIC)]
  if (verbose==T) {
    cat("  fit: lambda1",cs@par$lambda1[1],"lambda2",cs@par$lambda2[1],"\n")
    cat("  BIC = ",cs@par$value, "\n")
  }
  return(cs)
}

#' Check whether a normalization has converged
#' @export
#' 
has_converged = function(cs) {
  #return FALSE if legs have changed, and require at least 2 steps
  params=cs@diagnostics$params
  laststep=params[,step[.N]]
  if (laststep<=2) return(FALSE)
  if (!setequal(params[step==laststep,.(leg)],params[step==laststep-1,.(leg)])) return(FALSE)
  #check all legs present in the last step
  #check if all quantities directly involved in building the expected matrix have converged
  rel.precision=function(name,fn=identity){merge(params[step==laststep,.(leg,get(name))],
                               params[step==laststep-1,.(leg,get(name))],by="leg")[
                                 leg==leg[.N], max( abs(fn(V2.x[[1]])-fn(V2.y[[1]])) ) / ( max(fn(V2.x[[1]]))-min(fn(V2.x[[1]])) ) ]}
  conv.log_iota = rel.precision("log_iota")
  conv.log_rho = rel.precision("log_rho")
  conv.log_decay = rel.precision("log_decay")
  conv.signal = max(rel.precision("beta.phi"),rel.precision("lambda1"))
  #cat(" conv.log_iota ", conv.log_iota,
  #    " conv.log_rho ", conv.log_rho, " conv.log_decay ", conv.log_decay, " conv.phi ", conv.phi, "\n")
  cat(" relative precision for this iteration: iota ", conv.log_iota,
      " rho ", conv.log_rho, " decay ", conv.log_decay, " signal ", conv.signal, "\n")
  conv.param = all(c(conv.log_iota,conv.log_rho,
                     conv.log_decay,conv.signal)<cs@settings$tol, na.rm=T)
  return(conv.param)
}

#' count number of zeros in a given cut site, distance bin and signal bin
#' 
#' @return a data table, keyed on name,dbin,id1,pos1,bin1,bin2,dir,cat with the following columns
#' - name: name of the dataset
#' - id1, pos1, bin1: coordinates of the cut site
#' - bin2, dbin: signal/distance bin in which we looked at the intersections
#' - dir: fwd (rev) for contacts with downstream (resp upstream) cut-sites
#' - cat: whether we consider contacts on the left (contact L) or on the right (contact R) of this cut site
#' - ncross: number of cut site intersections (crossings) in this signal/distance/direction/category bin.
#'  We discard anything below cs@settings$dmin. Note that sum(ncross) is four times the total number of crossings (per dataset).
#' - nnz: number of non-zero contacts in this signal/distance/direction bin (max 2 per crossing).
#'  Note that sum(nnz) is twice the number of nonzeros (per dataset)
#' - nzero: number of zeros in this signal/distance/direction bin. We have nzero = 2*ncross - nnz.
#'  Note that sum(nzero) is twice the number of zeros (per dataset)
#'  For speed purposes downstream, we only return the entries where nzero>0, which corresponds to most of the entries anyway.
#'  We therefore have sum(nnz+nzero) equal to approximately twice the number of detectable counts, e.g. 8x the number of crossings
#'  
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
  chunksize=cs@biases[,ceiling(.N/(10*ncores))]
  nchunks=cs@biases[,ceiling(.N/chunksize)]
  registerDoParallel(cores=ncores)
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
  stopImplicitCluster()
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
  zeros=zeros[nzero>0]
  setkey(zeros,name,dbin,id1,pos1,bin1,bin2,dir,cat)
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
prepare_first_signal_estimation = function(biases, names, base.res) {
  ### build matrix
  #create an empty matrix containing all cells, even those with no cut-site intersection
  sbins=seq(biases[,min(pos)-1],biases[,max(pos)+1+base.res],base.res)
  signal.bins=unique(cut(c(sbins,head(sbins,n=-1)+base.res/2), sbins,
                         ordered_result=T, right=F, include.lowest=T,dig.lab=12))
  signal.mat=CJ(name=names,bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
  signal.mat[,c("phi","beta"):=list(0,0)]
  setkey(signal.mat,name,bin1,bin2)
  stopifnot(all(signal.mat[,.N,by=name]$N==signal.mat[,nlevels(bin1)*(nlevels(bin1)+1)/2]))
  return(list(signal=signal.mat,sbins=sbins))
}

#' Diagnostics plots to monitor convergence of normalization (gaussian
#' approximation)
#' 
#' @param cs a normalized CSnorm object
#' @param start,end the starting and ending steps on which to plot the info. Default is full range.
#'   
#' @return Two plots in a list. The first is for the four log-likelihoods, the
#'   second is for the parameters.
#' @export
#' 
#' @examples
plot_diagnostics = function(cs, start=1, end=Inf) {
  if (is.infinite(end) || cs@diagnostics$params[,max(step)]<end) end = cs@diagnostics$params[,max(step)] 
  plot=ggplot(cs@diagnostics$params[step>=start&step<end+1,.(step,leg,value)])+
    geom_line(aes(step,value))+geom_point(aes(step,value))+facet_wrap(~leg, scales = "free_y")+
    theme(legend.position="bottom")+scale_x_continuous(breaks = seq(start,end,2), minor_breaks = seq(start+1,end,2))
  vars=foreach(var=c("eC","eRJ","eDE","alpha","lambda_iota","lambda_rho","lambda_diag", "lambda1", "lambda2", "eCprime"),
               trans=(c(NA,NA,NA,NA,"log10","log10","log10","log10","log10",NA)),.combine=rbind) %do% get_all_values(cs,var,trans)
  plot2=ggplot(vars[step>=start&step<end+1])+geom_line(aes(step,value,group=as.integer(step)))+
    geom_point(aes(step,value,colour=leg))+facet_wrap(~variable, scales = "free_y")+
    scale_x_continuous(breaks = seq(start,end,2), minor_breaks = seq(start+1,end,2))
  return(list(plot=plot,plot2=plot2))
}

#' Save normalized CS object, discarding everything but
#' input counts and final normalization parameters
#' 
#' @param cs a normalized CSnorm object
#' @param fname where to save it to
#'   
#' @export
#' 
#' @examples
save_stripped = function(cs, fname) {
  cs@diagnostics=list(params=data.table())
  cs@zeros=data.table()
  strip_interaction = function(csi) {
    if (class(csi) %in% c("CSbsig","CSbdiff")) {
      csi@state=list()
      csi@cts=data.table()
    }
    return(csi)
  }
  strip_group = function(csg) {
    csg@cts=data.table()
    csg@interactions=lapply(csg@interactions,strip_interaction)
    return(csg)
  }
  cs@groups=lapply(cs@groups,strip_group)
  save(cs,file=fname)
}

#' Load normalized CS object saved with save_stripped,
#' rebuild zeros data.table
#' 
#' @param fname file name
#' @param ncores number of CPUs to use to reconstruct zeros
#'   
#' @return A cs object on which signal detection can be performed
#'
#' @export
#' 
#' @examples
load_stripped = function(fname, ncores=1) {
  cs=get(load(fname))
  cs@zeros = binless:::get_nzeros(cs, cs@settings$sbins, ncores=ncores)
  return(cs)
}


#' Compute initial exposures assuming a poisson model
#' @keywords internal
#' 
initial_guess_exposures = function(cs, cts.common, pseudocount=1e-2) {
  #biases
  cs@par$eDE=as.array(cs@biases[,log(pseudocount+mean(dangling.L+dangling.R)/2),keyby=c("name")]$V1)
  cs@par$eRJ=as.array(cs@biases[,log(pseudocount+mean(rejoined)),keyby=c("name")]$V1)
  cs@par$eC=array(0,cs@experiments[,.N])
  cs@par$eC = as.array(cts.common[,log(pseudocount+weighted.mean(count,weight)),keyby=c("name")]$V1)
  return(cs)
}

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

#' Cleanup a CSnorm object, store settings and populate it with initial guesses of all required parameters
#' @keywords internal
#' 
fresh_start = function(cs, bf_per_kb=50, bf_per_decade=10, bins_per_bf=10, base.res=5000,
                       bg.steps=5, iter=100, fit.signal=T, verbose=T, ncounts=100000, init.dispersion=1,
                       tol=1e-1, ncores=1, fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA) {
    #fresh start
    cs@par=list() #in case we have a weird object
    cs@groups=list()
    cs@diagnostics=list()
    #add settings
    if (fit.signal==F) base.res = cs@biases[,max(pos)-min(pos)]+2
    cs@settings = c(cs@settings[c("circularize","dmin","dmax","qmin","dfuse")],
                    list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf, base.res=base.res,
                         bg.steps=bg.steps, iter=iter, init.dispersion=init.dispersion, tol=tol,
                         fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                         fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at))
    cs@settings$Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    cs@settings$Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
    stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
    cs@settings$dbins=10**seq(log10(cs@settings$dmin-1),log10(cs@settings$dmax+1)+stepsz,stepsz)
    #initial guess
    if (verbose==T) cat("No initial guess provided\n")
    decay=CJ(name=cs@experiments[,name],dist=head(cs@settings$dbins,n=length(cs@settings$dbins)-1)*10**(stepsz/2))
    decay[,dbin:=cut(dist, cs@settings$dbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    decay[,c("dist","log_decay"):=list(NULL,0)]
    cs@par=list(eC=array(0,cs@experiments[,.N]), eRJ=array(0,cs@experiments[,.N]), eDE=array(0,cs@experiments[,.N]), alpha=init.dispersion,
                log_iota=array(0,cs@biases[,.N]), log_rho=array(0,cs@biases[,.N]),
                decay=decay, log_decay=0, tol_genomic=.1, tol_decay=.1, tol_disp=.1, tol_signal=1)
    #prepare signal matrix
    if (fit.signal==T) {
      if(verbose==T) cat("Preparing for signal estimation\n")
      stuff = binless:::prepare_first_signal_estimation(cs@biases, cs@experiments[,name], base.res)
      cs@par$signal=stuff$signal
      cs@par$beta.phi=stuff$signal[,beta]
      cs@par$lambda1=array(dim=cs@experiments[,.N])
      cs@settings$sbins=stuff$sbins
    } else {
      cs@settings$sbins=cs@biases[,c(min(pos)-1,max(pos)+1)]
      cs@par$signal=cs@biases[,.(phi=0,beta=0,bin1=cut(pos[1], cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                                 bin2=cut(pos[1], cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)),by=name]
      cs@par$beta.phi=cs@par$signal[,beta]
      cs@par$lambda1=NA
      setkey(cs@par$signal,name,bin1,bin2)
    }
    #get number of zeros along cut sites and decay
    if(verbose==T) cat("Counting zeros\n")
    cs@zeros = binless:::get_nzeros(cs, cs@settings$sbins, ncores=ncores)
    #set initial guess for exposures, decay and biases
    if(verbose==T) cat("Initial guess: residuals\n")
    cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
    if(verbose==T) cat("Initial guess: exposures\n")
    cs = binless:::initial_guess_exposures(cs, cts.common)
    if(verbose==T) cat("Initial guess: decay\n")
    cs = binless:::initial_guess_decay(cs, cts.common)
    if(verbose==T) cat("Initial guess: biases\n")
    cs = binless:::initial_guess_genomic(cs, cts.common)
    return(cs)
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
#' @param fit.signal boolean. Set to FALSE only for diagnostics.
#' @param verbose Display progress if TRUE
#' @param init.dispersion positive numeric. Value of the dispersion to use initially.
#' @param tol positive numeric (default 1e-2). Convergence tolerance on relative changes in the computed biases.
#' @param ncores positive integer (default 1). Number of cores to use.
#' @param fix.lambda1 whether to set lambda1 to a given value, or to estimate it
#' @param fix.lambda1.at if fix.lambda1==T, the approximate value where it is meant to be fixed. Might move a bit because
#'   of the positivity and degeneracy constraints.
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
#' 
normalize_binless = function(cs, restart=F, bf_per_kb=50, bf_per_decade=10, bins_per_bf=10, base.res=5000,
                     ngibbs = 15, bg.steps=5, iter=100, fit.signal=T,
                     verbose=T, ncounts=100000, init.dispersion=1,
                     tol=1e-1, ncores=1, fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #
  if (verbose==T) cat("Normalization with fast approximation and performance iteration\n")
  setkey(cs@biases, id, name)
  setkey(cs@counts, id1, id2, name)
  if (restart==F) {
    #fresh start
    cs = fresh_start(cs, bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, base.res = base.res,
                     bg.steps = bg.steps, iter = iter, fit.signal = fit.signal,
                     verbose = verbose, ncounts = ncounts, init.dispersion = init.dispersion,
                     tol = tol, ncores = ncores, fix.lambda1 = fix.lambda1, fix.lambda1.at = fix.lambda1.at,
                     fix.lambda2 = fix.lambda2, fix.lambda2.at = fix.lambda2.at)
    laststep=0
    #update eC everywhere in the first step
    update.eC=T
  } else {
    if (verbose==T) cat("Continuing already started normalization with its original settings\n")
    laststep = cs@diagnostics$params[,max(step)]
    update.eC = !(fit.signal==T && laststep > cs@settings$bg.steps)
  }
  #
  if(verbose==T) cat("Subsampling counts for dispersion\n")
  subcounts = binless:::subsample_counts(cs, ncounts)
  subcounts.weight = merge(cs@zeros[,.(nc=sum(ncross)/4),by=name],subcounts[,.(ns=.N),keyby=name],by="name")[,.(name,wt=nc/ns)]
  #gibbs sampling
  if (ngibbs==0) return(cs)
  for (i in (laststep + 1:ngibbs)) {
    if (verbose==T) cat("\n### Iteration",i,"\n")
    #
    #compute residuals once for this round
    if(verbose==T) cat(" Residuals\n")
    cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
    cs@diagnostics$residuals = rbind(cs@diagnostics$residuals,
                                     cts.common[bin1==min(bin1),.(step=i,z=weighted.mean(z,weight/var),phi=weighted.mean(phi,weight/var),
                                                                  log_decay=weighted.mean(log_decay,weight/var),log_bias=weighted.mean(log_bias,weight/var),
                                                                  log_mean=weighted.mean(lmu.nosig+phi,weight/var),weight=sum(weight/2),count=sum(weight*count/2)),
                                                keyby=c("name","bin2")])
    #fit iota and rho
    constrain.bias = fit.signal==T && i <= cs@settings$bg.steps+1
    a=system.time(cs <- binless:::gauss_genomic(cs, cts.common, update.eC=update.eC, verbose=verbose, constrain=constrain.bias))
    cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="bias", runtime=a[1]+a[4])
    #
    if (fit.signal==T && i > cs@settings$bg.steps) {
      if(verbose==T) cat(" Residuals\n")
      cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
      #fit signal using sparse fused lasso
      update.eC=F
      a=system.time(cs <- binless:::gauss_signal(cs, cts.common, verbose=verbose, ncores=ncores,
                                                       fix.lambda1=cs@settings$fix.lambda1,
                                                       fix.lambda1.at=cs@settings$fix.lambda1.at,
                                                       fix.lambda2=cs@settings$fix.lambda2,
                                                       fix.lambda2.at=cs@settings$fix.lambda2.at))
      cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="signal", runtime=a[1]+a[4])
    } else {
      #fit diagonal decay
      a=system.time(cs <- binless:::gauss_decay(cs, cts.common, update.eC=update.eC, verbose=verbose))
      cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="decay", runtime=a[1]+a[4])
    }
    #
    #fit dispersion
    a=system.time(cs <- binless:::gauss_dispersion(cs, counts=subcounts, weight=subcounts.weight, verbose=verbose))
    cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="disp", runtime=a[1]+a[4])
    #
    #check for convergence
    if (has_converged(cs)) {
      if (fit.signal == T && i <= cs@settings$bg.steps) {
        cat("Background has converged, fitting signal\n")
        cs@settings$bg.steps = i #compute signal at next step
      } else {
        if (verbose==T) {
          cat("Normalization has converged\n")
          cat(" setting parameters to their unconstrained values\n")
        }
        if (fit.signal==T) {
          setnames(cs@par$signal,c("phi","phi.unconstr"),c("phi.constr","phi"))
          cs@par$eC = cs@par$eC - cs@par$eCprime + cs@par$eCprime.unconstr
        }
        break
      }
    }
  }
  if (verbose==T) cat("Done\n")
  return(cs)
}


