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
                 ngibbs = 25, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
                 ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
} else {
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,".RData"))
  cs = run_gauss(cs, restart=T, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
                 ngibbs = 10, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
                 ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
}
save(cs,file=paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,".RData"))

if (F) {
  
  base.res=5000
  load(paste0("data/rao_HiCall_",sub,"_IMR90_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,".RData"))
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
  ggplot(signals[step>=step[.N]-1])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_grid(step~ name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phihat))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  
  subs=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "Fig1C_1M","FOXP1_1.3M","TBX3_1.5M",
         "Comparison_1.7M","22qter_1.7M", "Talk_2M", "ADAMTS1_2.3M")
  info = foreach(sub=subs,.combine=rbind) %dopar% {
    #info = foreach(sub=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "Fig1C_1M"),.combine=rbind) %dopar% {
    load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk30_dfuse5qmin_0.05_cv_cvsd_outlier_rmdiag_simplified_bounds.RData"))
    #load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"qmin_",qmin,"_cv_cvsd_outlier_rmdiag_simplified.RData"))
    #load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"_",run,".RData"))
    data.table(ori=sub,has.converged=csnorm:::has_converged(cs),#run=run,
               runtime=cs@diagnostics$params[,sum(runtime)],nsteps=cs@diagnostics$params[,max(step)],
               name=cs@experiments[,name],
               lambda1=cs@par$lambda1,lambda2=cs@par$lambda2,eCprime=cs@par$eCprime,lambda_diag=cs@par$lambda_diag)#,
    #lambda1.post=cs@groups[[1]]@interactions[[1]]@par$lambda1, lambda2.post=cs@groups[[1]]@interactions[[1]]@par$lambda2)
  }
  info
  
  
  rtinfo = foreach(sub=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "Fig1C_1M","FOXP1_1.3M","TBX3_1.5M",
                         "Comparison_1.7M","22qter_1.7M", "Talk_2M", "ADAMTS1_2.3M"),.combine=rbind) %dopar% {
                           load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk30_dfuse5qmin_0.05_cv_cvsd_outlier_rmdiag_simplified.RData"))
                           cs@diagnostics$params[,.(sub,step,leg,runtime)]
                         }
  ggplot(rtinfo)+geom_line(aes(step,runtime,colour=leg))+scale_y_log10()+facet_wrap(~sub)
  
}

