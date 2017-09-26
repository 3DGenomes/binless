library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

args=commandArgs(trailingOnly=TRUE)
sub="Tbx19_700k"
bpk=30
dfuse=5 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
ncores=4
qmin=0.05
base.res=as.integer(args[1])

setwd("/home/yannick/simulations/cs_norm")

restart=F
if (restart==F) {
  load(paste0("data/rao_HiCall_IMR90_",sub,"_csdata.RData"))
  csd2=csd
  cs=merge_cs_norm_datasets(list(csd2), different.decays="none", dfuse=dfuse, qmin=qmin)
  cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
                 ngibbs = 20, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
                 ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
} else {
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_grpall.RData"))
  cs = run_gauss(cs, restart=T, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
                 ngibbs = 10, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
                 ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
}
save(cs,file=paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_grpall.RData"))

if (F) {
  
  base.res=2500
  #load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,".RData"))
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_grpall.RData"))
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_newdisp.RData"))
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_newdisp_grpall.RData"))
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_newdisp_grpall_lindecay.RData"))
  csnorm:::has_converged(cs)
  cs@diagnostics$params[,sum(runtime)/3600]
  ggplot(cs@diagnostics$params[,.(step,leg,runtime)])+geom_line(aes(step,runtime,colour=leg))+scale_y_log10()
  
  plot_diagnostics(cs)$plot
  plot_diagnostics(cs)$plot2
  
  signals=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
    if ("signal" %in% cs@diagnostics$params[step==i,leg]) {
      sig=copy(cs@diagnostics$params[step==i&leg=="signal",signal][[1]])
      sig[,step:=i]
      sig
    }
  }
  
  ggplot(signals[name==name[1]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[name==name[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  #ggplot(signals[name==name[1]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(data=signals.old[name==name[1]],aes(bin2,bin1,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]-1])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_grid(step~ name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=pmax(-1,pmin(phihat,max(phi)))))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  
  info = foreach(base.res=c(10000,5000,2500),.combine=rbind, .errorhandling="remove") %:%
         foreach(ext=c("","_grpall","_newdisp","_newdisp_grpall","_newdisp_grpall_lindecay"),.combine=rbind) %dopar% {
    load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,ext,".RData"))
    data.table(ori=sub,base.res=base.res,ext=ext,runtime=cs@diagnostics$params[,sum(runtime)],
               nsteps=cs@diagnostics$params[,max(step)],name=cs@experiments[,name],
               lambda1=cs@par$lambda1,lambda2=cs@par$lambda2,eCprime=cs@par$eCprime,lambda_diag=cs@par$lambda_diag)
  }
  info
  
  mats = foreach(base.res=c(10000,5000,2500),.combine=rbind, .errorhandling="remove") %:%
    foreach(ext=c("","_grpall","_newdisp","_newdisp_grpall","_newdisp_grpall_lindecay"),.combine=rbind) %dopar% {
      load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,ext,".RData"))
      mat=cs@par$signal
      bin1.begin=mat[,bin1]
      bin2.begin=mat[,bin2]
      levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
      levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
      mat[,begin1:=as.integer(as.character(bin1.begin))]
      mat[,begin2:=as.integer(as.character(bin2.begin))]
      #ggplot(mat)+geom_raster(aes(begin2,begin1,fill=phihat))+geom_raster(aes(begin1,begin2,fill=phi))+scale_fill_gradient2()+facet_wrap(~name)+coord_fixed()
      mat[,c("ori","base.res","ext"):=list(sub,base.res,ext)]
    }
  ggplot(mats)+geom_raster(aes(begin2,begin1,fill=pmin(phihat,3)))+geom_raster(aes(begin1,begin2,fill=phi))+
    scale_fill_gradient2()+facet_grid(ext~base.res)+coord_fixed()
  ggplot(mats)+geom_raster(aes(begin2,begin1,fill=pmin(phi,0.5)))+geom_raster(aes(begin1,begin2,fill=pmin(phi,0.5)))+
    scale_fill_gradient2()+facet_grid(ext~base.res)+coord_fixed()
  
}

