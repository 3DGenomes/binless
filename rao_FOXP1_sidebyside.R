library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

bpk=30
dfuse=20 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
ncores=5
qmin=0.01
base.res=10000

setwd("/home/yannick/simulations/cs_norm")

sub="FOXP1_1.3M"
load(paste0("data/rao_HiCall_GM12878_",sub,"_csdata.RData"))
csd1=csd
load(paste0("data/rao_HiCall_IMR90_",sub,"_csdata.RData"))
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none", dfuse=dfuse, qmin=qmin)
cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 25, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
               ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag.RData"))

for (resolution in c(5000,10000,20000,50000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_binned_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
  cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
  save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag.RData"))
}

if (F) {
  
  #check all runs, pick best
  load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_base10k_bpk30_dfuse20_cv_cvsd_outlier_rmdiag.RData")
  cs10=cs
  load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_base5k_bpk50_dfuse5qmin_0.05_cv_cvsd_outlier_rmdiag.RData")
  cs5=cs
  load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_base20k_bpk50_dfuse5qmin_0.05_cv_cvsd_outlier_rmdiag.RData")
  cs20=cs
  load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_base5k_bpk30_dfuse5_cv_cvsd_outlier_rmdiag.RData")
  cs5b=cs
  
  #info
  info = foreach(run=c("cs5","cs5b","cs10","cs20"),cs=c(cs5,cs5b,cs10,cs20), .combine=rbind) %do% {
    data.table(run=run,qmin=cs@settings$qmin,dfuse=cs@settings$dfuse,base.res=cs@settings$base.res,
               cv=csnorm:::has_converged(cs), resolutions=list(sapply(cs@groups,function(x){x@resolution})))
  }
  
  #convergence
  mat = foreach(run=c("cs5","cs5b","cs10","cs20"),cs=c(cs5,cs5b,cs10,cs20), .combine=rbind) %do% {
    sig=rbindlist(cs@diagnostics$params[leg=="signal",signal],use.names = T, idcol = "step")[step>=max(step)-1]
    sig[,run:=run]
  }
  ggplot(mat[name==name[.N]])+
    geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()+
    geom_raster(aes(bin2,bin1,fill=phi))+coord_fixed()+facet_grid(run~name+step)
  
  #binless signal
  mat = foreach(run=c("cs5","cs5b","cs10","cs20"),cs=c(cs5,cs5b,cs10,cs20), .combine=rbind) %do% {
    sig=get_interactions(cs, resolution=5000, group="all",type="CSbsig")
    sig[,run:=run]
  }
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=phi))+scale_fill_gradient2()+
    geom_raster(aes(begin2,begin1,fill=phi))+coord_fixed()+facet_grid(run~name)
  
  #binless differences
  mat = foreach(run=c("cs5","cs5b","cs10","cs20"),cs=c(cs5,cs5b,cs10,cs20), .combine=rbind) %do% {
    sig=get_interactions(cs, resolution=5000, group="all",type="CSbdiff", ref=cs@experiments[1,name])
    sig[,run:=run]
  }
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=delta))+scale_fill_gradient2()+
    geom_raster(aes(begin2,begin1,fill=delta))+coord_fixed()+facet_grid(run~name)
  
  
  
  #pick dataset to show 4 different resolutions
  load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_base10k_bpk30_dfuse20_cv_cvsd_outlier_rmdiag.RData")
  
  mat = foreach(resolution=c(5000,10000,20000,50000), .combine=rbind, .errorhandling="remove") %do% {
    mat=get_matrices(cs, resolution=resolution, group="all")[name==name[1]]
    mat[,resolution:=resolution]
  }
  melted.mat=melt(mat[,.(resolution,name,begin1,begin2,observed,decay=decaymat,genomic=biasmat,
                         normalized,normalized.sd,
                         signal,signal.sd)],
                  id.vars=c("name","begin1","begin2","resolution"))
  #melted.mat[,value:=value/.SD[begin1==begin2,max(value)],by=c("variable","name")]
  melted.mat[,value:=pmin(value,quantile(value,0.99,na.rm=T)),by=c("variable","name","resolution")]
  melted.mat[,value:=(value-min(value,na.rm=T))/(max(value,na.rm=T)-min(value,na.rm=T)),by=c("variable","name","resolution")]
  #one dataset, all bw matrices
  ggplot(melted.mat)+geom_raster(aes(begin1,begin2,fill=value))+
    geom_raster(aes(begin2,begin1,fill=value))+facet_grid(resolution~variable)+scale_fill_gradient(low="white",high="black",na.value="white")+
    coord_fixed()+theme(panel.background=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="genomic coordinate",y="genomic coordinate")+guides(fill=F)
  ggsave(filename=paste0("images/rao_FOXP1_4resolutions_base10_binned.png"),width=32,height=20)
  

  #side-by-side at 5k: observed
  resolution=5000
  mat=get_matrices(cs, resolution=resolution, group="all")
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=observed))+
    geom_raster(aes(begin2,begin1,fill=observed))+coord_fixed()+facet_wrap(~name)+
    scale_fill_gradient(low="white",high="black",na.value="white",trans="log")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/rao_FOXP1_GMvsIMR90_base10at5_observed.png"),width=10,height=5)
  
  #side-by-side at 5k: signal
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSsig")
  #mat[,is.significant:=prob.gt.expected>1-1e-6]
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=log10(signal)))+
    geom_raster(aes(begin2,begin1,fill=log10(signal)))+coord_fixed()+facet_wrap(~name)+
    geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm")) + labs(fill="log10 FC")
  ggsave(filename=paste0("images/rao_FOXP1_GMvsIMR90_base10at5_binned_signal.png"),width=10,height=5)
  
  #side-by-side at 5k: differences
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSdiff")
  mat[,name:="FOXP1 IMR90 vs GM"]
  #mat[,is.significant:=prob.gt.expected>1-1e-6]
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=log10(difference)))+
    geom_raster(aes(begin2,begin1,fill=log10(difference)))+coord_fixed()+facet_wrap(~name)+
    geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm")) + labs(fill="log10 FC")
  ggsave(filename=paste0("images/rao_FOXP1_GMvsIMR90_base10at5_binned_difference.png"),width=8,height=5)
  
}
