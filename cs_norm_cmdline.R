library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)
library(scam)
library(Hmisc)
library(foreach)
library(doParallel)
registerDoParallel(cores=30)

setwd("/home/yannick/simulations/cs_norm")


fill_zeros = function(biases,counts) {
  newcounts=CJ(biases[,paste(id,pos)],biases[,paste(id,pos)])
  newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
  newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
  newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
  newcounts=newcounts[pos1<pos2]
  setkey(newcounts, id1, id2, pos1, pos2)
  setkey(counts, id1, id2, pos1, pos2)
  newcounts=counts[newcounts]
  newcounts[is.na(contact.close),contact.close:=0]
  newcounts[is.na(contact.far),contact.far:=0]
  newcounts[is.na(contact.up),contact.up:=0]
  newcounts[is.na(contact.down),contact.down:=0]
  return(newcounts)
}

optimize_all_meanfield = function(model, biases, counts, meanfield, maxcount, bf_per_kb=1, bf_per_decade=5,
                                  iter=10000, verbose=T, init=0, ...) {
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  cclose=counts[contact.close>maxcount,.(id1,id2,distance,count=contact.close)]
  cfar=counts[contact.far>maxcount,.(id1,id2,distance,count=contact.far)]
  cup=counts[contact.up>maxcount,.(id1,id2,distance,count=contact.up)]
  cdown=counts[contact.down>maxcount,.(id1,id2,distance,count=contact.down)]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               dmin=dmin, dmax=dmax, Kdiag=Kdiag,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               Nl=mf$Nkl[,.N], Nkl_count=mf$Nkl[,count], Nkl_cidx=mf$Nkl[,id], Nkl_N=mf$Nkl[,N], Nkl_levels=mf$Nkl[,sum(diff(N)!=0)+1],
               Nr=mf$Nkr[,.N], Nkr_count=mf$Nkr[,count], Nkr_cidx=mf$Nkr[,id], Nkr_N=mf$Nkr[,N], Nkr_levels=mf$Nkr[,sum(diff(N)!=0)+1],
               Nd=mf$Nkd[,.N], Nkd_count=mf$Nkd[,count], Nkd_d=mf$Nkd[,mdist], Nkd_N=mf$Nkd[,N], Nkd_levels=mf$Nkd[,sum(diff(N)!=0)+1])
  if (data$Nl==0) data$Nkl_levels=0
  if (data$Nr==0) data$Nkr_levels=0
  if (data$Nd==0) data$Nkd_levels=0
  message("Mean field optimization")
  message("Krow        : ", Krow)
  message("Kdiag       : ", Kdiag)
  message("Biases      : ", biases[,.N])
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  message("Left counts : ", mf$Nkl[,.N], " (", data$Nkl_levels, " levels)")
  message("Right counts: ", mf$Nkr[,.N], " (", data$Nkr_levels, " levels)")
  message("Decay counts: ", mf$Nkd[,.N], " (", data$Nkd_levels, " levels)")
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

predict_all_meanfield = function(model, biases, counts, meanfield, opt, bf_per_decade=5, verbose=T) {
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag_centered=opt$par$beta_diag_centered)
  message("Mean field prediction")
  message("Kdiag       : ", Kdiag)
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
}

get_binned_matrices = function(model, biases, counts, meanfield, opt, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL, bf_per_decade=5, verbose=T) {
  stopifnot(counts[distance!=pos2-pos1,.N]==0) #need to implement circular genomes
  csub=copy(counts) #need to implement taking only needed part of matrix, and reporting log_nu and log_delta appropriately
  bsub=copy(biases)
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  Kdiag=round(csub[,(log10(dmax)-log10(dmin))*bf_per_decade])
  npoints=100*Kdiag #evaluate spline with 100 equidistant points per basis function
  #bin existing counts and biases
  if (is.null(b1)) b1=bsub[,min(pos)]-1
  if (is.null(b2)) b2=bsub[,min(pos)]-1
  if (is.null(e1)) e1=bsub[,max(pos)]+1
  if (is.null(e2)) e2=bsub[,max(pos)]+1
  bins1=seq(b1,e1+resolution,resolution)
  bins2=seq(b2,e2+resolution,resolution)
  csub[,c("bin1","bin2"):=list(cut(pos1, bins1, ordered_result=T, right=F, include.lowest=T),
                                 cut(pos2, bins2, ordered_result=T, right=F, include.lowest=T))]
  csub[,log_mean_cup:=opt$par$log_mean_cup]
  csub[,log_mean_cdown:=opt$par$log_mean_cdown]
  csub[,log_mean_cclose:=opt$par$log_mean_cclose]
  csub[,log_mean_cfar:=opt$par$log_mean_cfar]
  bsub[,log_nu:=opt$par$log_nu]
  bsub[,log_delta:=opt$par$log_delta]
  bsub[,c("bin1","bin2"):=list(cut(pos, bins1, ordered_result=T, right=F, include.lowest=T),
                                 cut(pos, bins2, ordered_result=T, right=F, include.lowest=T))]
  bsub[,c("ibin1","ibin2"):=list(as.integer(bin1),as.integer(bin2))]
  #run computation of mean
  cclose=csub[,.(id1,id2,bin1,bin2,distance,count=contact.close,logmean=log_mean_cclose)]
  cfar=csub[,.(id1,id2,bin1,bin2,distance,count=contact.far,logmean=log_mean_cfar)]
  cup=csub[,.(id1,id2,bin1,bin2,distance,count=contact.up,logmean=log_mean_cup)]
  cdown=csub[,.(id1,id2,bin1,bin2,distance,count=contact.down,logmean=log_mean_cdown)]
  csub=rbind(cclose,cfar,cup,cdown)[!(is.na(bin1)|is.na(bin2)|count==0)]
  data = list( Kdiag=Kdiag, 
               S1=bsub[!is.na(bin1),.N], S2=bsub[!is.na(bin2),.N], 
               cutsites1=bsub[!is.na(bin1),pos], cutsites2=bsub[!is.na(bin2),pos],
               dmin=dmin, dmax=dmax, npoints=npoints,
               N=csub[,.N],   counts=csub[,count], cdist=csub[,distance], cmean=csub[,exp(logmean)],
               eC=opt$par$eC,
               log_nu1=bsub[!is.na(bin1),log_nu], log_nu2=bsub[!is.na(bin2),log_nu],
               log_delta1=bsub[!is.na(bin1),log_delta], log_delta2=bsub[!is.na(bin2),log_delta],
               beta_diag_centered=opt$par$beta_diag_centered,
               B1=csub[,nlevels(bin1)], B2=csub[,nlevels(bin2)],
               cbins1=csub[,as.integer(bin1)], cbins2=csub[,as.integer(bin2)],
               bbins1=bsub[!is.na(bin1),ibin1], bbins2=bsub[!is.na(bin2),ibin2])
  binned=optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
  #format output
  so = data.table(melt(binned$par$observed))
  setnames(so,"value","observed")
  so[,expected:=melt(binned$par$expected)$value]
  so[,ncounts:=melt(binned$par$ncounts)$value]
  so[,normalized:=melt(binned$par$normalized)$value]
  so=so[ncounts>0]
  setkey(so,Var1)
  so=bsub[,.(bin1=bin1[1]),keyby=ibin1][so]
  setkey(so,Var2)
  so=bsub[,.(bin2=bin2[1]),keyby=ibin2][so]
  so[,c("ibin1","ibin2"):=list(NULL,NULL)]
  so[,lFC:=log2(observed/expected)]
  return(list(distance=exp(binned$par$log_dist), decay=exp(binned$par$log_decay), mat=so))
}

get_dispersions = function(model, binned, iter=10000, verbose=T) {
  data=list(B=binned[,.N],observed=binned[,observed],expected=binned[,expected],ncounts=binned[,ncounts])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
}

#compute p(gamma2>gamma1) = \int_{0}^{+infty} dx p_gamma2(x) \int_{0}^{x} dy p_gamma1(y)
compute_gamma_overlap = function(alpha1,beta1,alpha2,beta2, bounds=5, ncores=1) {
  registerDoParallel(cores=ncores)
  foreach (a1=alpha1, b1=beta1, a2=alpha2, b2=beta2, .packages="stats", .combine=c) %dopar% {
    #infinite integral does not work very well so truncate outer integral
    mu1=a1/b1
    mu2=a2/b2
    sd1=sqrt(a1)/b1
    sd2=sqrt(a2)/b2
    xmin = max(0,min(mu1-bounds*sd1,mu2-bounds*sd2))
    xmax = max(mu1+bounds*sd1,mu2+bounds*sd2)
    a=integrate(function(x){exp(dgamma(x,a2,rate=b2,log=T)+pgamma(x,a1,rate=b1,log.p=T))},xmin,xmax)
    if (a$abs.error<=0) {NA} else {a$value}
  }
}

detect_interactions = function(binned, dispersions, threshold=0.95, ncores=1){
  #report gamma parameters
  mat=copy(binned)
  mat[,c("alpha1","beta1"):=list(dispersions,dispersions/expected)]
  mat[,c("alpha2","beta2"):=list(alpha1+observed,beta1+1)]
  mat[,prob.observed.gt.expected:=compute_gamma_overlap(alpha1,beta1,alpha2,beta2,ncores=ncores)]
  mat[,is.interaction:=prob.observed.gt.expected>threshold | 1-prob.observed.gt.expected>threshold]
  #write begins/ends
  bin1.begin=mat[,bin1]
  bin1.end=mat[,bin1]
  bin2.begin=mat[,bin2]
  bin2.end=mat[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[[,]")[2][[1]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[[,)]")[2][[1]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[[,]")[2][[1]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[[,)]")[2][[1]]
  mat[,begin1:=as.integer(as.character(bin1.begin))]
  mat[,end1:=as.integer(as.character(bin1.end))]
  mat[,begin2:=as.integer(as.character(bin2.begin))]
  mat[,end2:=as.integer(as.character(bin2.end))]
  return(mat)
}

#estimates the values of the count or dispersion required to cross a given threshold
thresholds_estimator = function(observed, expected, dispersion, threshold=0.95, counts.range=c(expected/100,expected*100), disp.range=c(0,100*dispersion), compute=T) {
  #plot gamma distributions
  a1=dispersion
  b1=dispersion/expected
  a2=a1+observed
  b2=b1+1
  mu1=a1/b1
  mu2=a2/b2
  sd1=sqrt(a1)/b1
  sd2=sqrt(a2)/b2
  message("prior:     mean ",mu1, " sd ",sd1)
  message("posterior: mean ",mu2, " sd ",sd2)
  xmin = min(c(0,counts.range))
  xmax = max(counts.range)
  message("evaluation bounds: [",xmin,",",xmax,"]")
  p1=ggplot(data.table(x=c(xmin,xmax)),aes(x))+
    stat_function(fun=function(x){dgamma(x,a1,rate=b1)},aes(colour="prior"))+
    stat_function(fun=function(x){dgamma(x,a2,rate=b2)},aes(colour="posterior"))+
    scale_colour_manual("Gamma", values = c("blue", "red"))
  #helper function  
  gammaQuery=function(mu,theta,count){
    alpha1=theta
    alpha2=theta+count
    beta1=theta/mu
    beta2=beta1+1
    return(function(x){exp(dgamma(x,alpha2,rate=beta2,log=T)+pgamma(x,alpha1,rate=beta1,log.p=T))})
  }
  #query counts
  igammaQuery_c=function(x){integrate(gammaQuery(expected,dispersion,count=x),xmin,xmax)$value-threshold}
  p2=ggplot(data.table(x=c(xmin,xmax),y=c(0,1)),aes(x))+
    stat_function(fun=Vectorize(function(x){igammaQuery_c(x)+threshold}))+
    geom_hline(aes(yintercept=threshold))+ geom_hline(aes(yintercept=1-threshold))+
    labs(x="observed count", y="P(observed>expected)")
  if (compute){
    clo = uniroot(function(x){igammaQuery_c(x)+2*threshold-1},c(xmin,expected))$root
    chi = uniroot(igammaQuery_c,c(expected,xmax))$root
    message("At a threshold of ",threshold," significant counts are lower than ",clo," and higher than ",chi)
    p2=p2+geom_vline(aes(xintercept=clo))+ geom_vline(aes(xintercept=chi))
  }
  #query dispersion
  igammaQuery_theta=function(x){integrate(gammaQuery(expected,x,observed),xmin,xmax)$value-threshold}
  p3=ggplot(data.table(x=disp.range,y=c(0,1)),aes(x))+
    stat_function(fun=Vectorize(function(x){igammaQuery_theta(x)+threshold}))+
    geom_hline(aes(yintercept=threshold))+geom_hline(aes(yintercept=1-threshold))+
    labs(x="dispersion", y="P(observed>expected)")+scale_x_log10()
  if (compute) {
    dstar = uniroot(igammaQuery_theta,disp.range)$root
    p3=p3+geom_vline(aes(xintercept=dstar))
    if (observed > expected) {
      message("At a threshold of ",threshold," the dispersion must be larger than ",dstar," to reach significance")
    } else {
      message("At a threshold of ",threshold," the dispersion must be smaller than ",dstar," to reach significance")
    }
  }
  return(list(p1,p2,p3))
}

read_and_prepare = function(infile, outprefix, skip=0L, both=T, distance_bins_per_decade=100, circularize=-1) {
  message("*** READ")
  data=read_tsv(infile, skip=skip)
  message("*** CATEGORIZE")
  data = categorize_by_new_type(data)
  message("*** BIASES AND COUNTS")
  cs_data = prepare_for_sparse_cs_norm(data, both=both, circularize=circularize)
  dset_statistics(cs_data$biases,cs_data$counts)
  message("*** MEANFIELD")
  cs_data$meanfield=bin_for_mean_field(cs_data$biases, cs_data$counts, distance_bins_per_decade=distance_bins_per_decade)
  message("*** WRITE")
  write.table(cs_data$biases, file = paste0(outprefix,"_biases.dat"), quote=F, row.names = F)
  write.table(cs_data$counts, file = paste0(outprefix,"_counts.dat"), quote=F, row.names = F)
  if (both==T) {
    both=cs_data$both
    save(both, file = paste0(outprefix,"_both.RData"))
  }
  meanfield=cs_data$meanfield
  save(meanfield, file = paste0(outprefix, "_meanfield_", distance_bins_per_decade, ".RData"))
  return(cs_data)
}




args=commandArgs(trailingOnly = T)
prefix=args[1]
if (length(args)>=2) suffix=args[2] else suffix=NULL

message("normalization on ",prefix)

biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))
#meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)
#save(meanfield, file = paste0("data/",prefix,"_meanfield_100.RData"))
load(paste0("data/",prefix,"_meanfield_100.RData"))

#### optimization wihout prior guesses
smfit = stan_model(file = "cs_norm_fit.stan")
smpred = stan_model(file = "cs_norm_predict.stan")
smbin = stan_model("cs_norm_predict_binned.stan")
smdisp = stan_model("cs_norm_binned_dispersions.stan")
maxcount=-1
a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                           bf_per_decade=5, verbose = T, iter=100000)) #, tol_rel_grad=1e3, tol_rel_obj=1e3
op$pred=predict_all_meanfield(smpred, biases, counts, meanfield, op, verbose=T)$par
op$binned=get_binned_matrices(smbin, biases, counts, meanfield, op, resolution=1000)
op$disp=get_dispersions(smdisp, op$binned$mat)$par
op$mat=detect_interactions(op$binned$mat, op$disp$dispersion, ncores=30) #interaction detection using binned dispersion estimates
save(op, file=paste0("data/",prefix,"_op_maxcount_",maxcount,".RData"))


