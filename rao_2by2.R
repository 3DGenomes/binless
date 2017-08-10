library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)

args=commandArgs(trailingOnly=TRUE)
sub=args[1]
bpk=30
dfuse=5 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
ncores=5
qmin=0.01
base.res=10000

setwd("/home/yannick/simulations/cs_norm")

load(paste0("data/rao_HiC003_GM12878_",sub,"_csdata.RData"))
csd1=csd
load(paste0("data/rao_HiC006_GM12878_",sub,"_csdata.RData"))
csd2=csd
load(paste0("data/rao_HiC053_IMR90_",sub,"_csdata.RData"))
csd3=csd
load(paste0("data/rao_HiC055_IMR90_",sub,"_csdata.RData"))
csd4=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3,csd4), different.decays="none", dfuse=dfuse, qmin=qmin)
cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 20, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
               ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k_dfuse",dfuse,"_qmin_",qmin,".RData"))


foreach (resolution=c(5000,10000,20000),.errorhandling = "pass") %do% {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_binned_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
  cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
  save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k_dfuse",dfuse,"_qmin_",qmin,".RData"))
}

foreach (resolution=c(5000,10000,20000),.errorhandling = "pass") %do% {
  cs=group_datasets(cs, resolution=resolution, verbose=T, ncores=ncores, group="condition")
  cs=detect_binned_interactions(cs, resolution=resolution, group="condition", ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="condition", ncores=ncores)
  ref=get_matrices(cs,resolution,group="condition")[,name[1]]
  cs=detect_binned_differences(cs, resolution=resolution, group="condition", ncores=ncores, ref=ref)
  cs=detect_binless_differences(cs, resolution=resolution, group="condition", ncores=ncores, ref=ref)
  save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k_dfuse",dfuse,"_qmin_",qmin,".RData"))
}

if (F) {
  subs=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "FOXP1_1.3M",
         "TBX3_1.5M","Comparison_1.7M")
  
  #plot last signal
  info = foreach(sub=subs, .combine=rbind, .errorhandling="remove") %dopar% {
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base10k_dfuse",5,"_qmin_",0.01,"_stripped.RData"))
    ggplot(cs@par$signal)+geom_raster(aes(bin1,bin2,fill=phi))+
      geom_raster(aes(bin2,bin1,fill=pmin(phihat,3)))+facet_wrap(~name)+
      scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
      theme_void()+ theme(axis.title=element_blank(),
                          panel.background = element_rect(fill = "white", colour = "black"),
                          panel.spacing=unit(0,"cm"))+
      coord_fixed()+labs(fill="log10 FC")
    ggsave(filename=paste0("images/rao_2by2_",sub,"_base10_iterated_signal.png"),width=10,height=10)
  }
  
  #get binless parameters
  info = foreach(sub=subs, .combine=rbind) %dopar% {
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base10k_dfuse",5,"_qmin_",0.01,"_stripped.RData"))
    foreach (resolution=c(5000,10000), .combine=rbind, .errorhandling = "remove") %do% {
      #side-by-side at 5k: observed
      idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
      csg=cs@groups[[idx1]]
      idx2=get_cs_interaction_idx(csg, type="CSbsig", raise=T)
      csi=csg@interactions[[idx2]]
      ret=as.data.table(csi@par)
      ret[,c("run","resolution"):=list(sub,resolution)]
      ret
    }
  }
  
  
  #generate plots
  mat = foreach(sub=subs, .combine=rbind, .errorhandling="remove") %dopar% {
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base10k_dfuse",5,"_qmin_",0.01,"_stripped.RData"))
    foreach (resolution=c(5000,10000), .errorhandling="remove") %do% {
      #side-by-side at 5k: observed
      mat=get_matrices(cs, resolution=resolution, group="all")
      ggplot(mat)+
        geom_raster(aes(begin1,begin2,fill=observed))+
        geom_raster(aes(begin2,begin1,fill=observed))+coord_fixed()+facet_wrap(~name)+
        scale_fill_gradient(low="white",high="black",na.value="white",trans="log")+
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm"))
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_observed.png"),width=10,height=10)
      
      #side-by-side at 5k: signal
      mat=get_interactions(cs, resolution=resolution, group="all", type="CSsig")
      ggplot(mat)+
        geom_raster(aes(begin1,begin2,fill=log10(signal)))+
        geom_raster(aes(begin2,begin1,fill=log10(signal)))+coord_fixed()+facet_wrap(~name)+
        geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+
        scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm")) + labs(fill="log10 FC")
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binned_signal.png"),width=10,height=10)
      
      #side-by-side at 5k: differences
      mat=get_interactions(cs, resolution=resolution, group="all", type="CSdiff")
      ggplot(mat)+
        geom_raster(aes(begin1,begin2,fill=log10(difference)))+
        geom_raster(aes(begin2,begin1,fill=log10(difference)))+coord_fixed()+facet_wrap(~name)+
        geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+
        scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm")) + labs(fill="log10 FC")
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binned_differences.png"),width=15,height=5)
      
      #side-by-side at 5k: binless signal without threshold
      mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
      ggplot(mat)+geom_raster(aes(begin1,begin2,fill=beta))+
        geom_raster(aes(begin2,begin1,fill=beta_cv))+facet_wrap(~name)+
        scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm"))+
        coord_fixed()+labs(fill="log10 FC")
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binless_signal_nothresh.png"),width=10,height=10)
      
      #side-by-side at 5k: binless signal
      idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
      csg=cs@groups[[idx1]]
      idx2=get_cs_interaction_idx(csg, type="CSbsig", raise=T)
      csi=csg@interactions[[idx2]]
      csi@settings$min.l10FC=0.2
      mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
      mat[,value:=phi]
      mat = csnorm:::detect_binless_patches(mat, csi@settings)
      mat[,value:=NULL]
      mat[,phi.max:=ifelse(is.maximum==T,NA,phi)]
      ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi))+
        geom_raster(aes(begin2,begin1,fill=phi.max))+facet_wrap(~name)+
        scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="black")+
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm"))+
        coord_fixed()+labs(fill="log10 FC")
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binless_signal.png"),width=10,height=10)
      
      
      #side-by-side at 5k: binless difference without threshold
      mat=get_interactions(cs, type="CSbdiff", resolution=resolution, group="all", ref=cs@experiments[1,name])
      ggplot(mat)+geom_raster(aes(begin1,begin2,fill=beta))+
        geom_raster(aes(begin2,begin1,fill=beta_cv))+facet_wrap(~name)+
        scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white") +
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm"))+
        coord_fixed()+labs(fill="log10 FC")
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binless_differences_nothresh.png"),width=15,height=5)
      
      #side-by-side at 5k: binless difference
      idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
      csg=cs@groups[[idx1]]
      idx2=get_cs_interaction_idx(csg, type="CSbdiff", raise=T, ref=cs@experiments[1,name])
      csi=csg@interactions[[idx2]]
      csi@settings$min.l10FC=0.2
      mat=get_interactions(cs, type="CSbdiff", resolution=resolution, group="all", ref=cs@experiments[1,name])
      mat[,value:=delta]
      mat = csnorm:::detect_binless_patches(mat, csi@settings)
      mat[,value:=NULL]
      mat[,delta.opt:=ifelse(is.maximum==T | is.minimum==T,NA,delta)]
      ggplot(mat)+geom_raster(aes(begin1,begin2,fill=delta))+
        geom_raster(aes(begin2,begin1,fill=delta.opt))+facet_wrap(~name)+
        scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="black") +
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
        theme_void()+ theme(axis.title=element_blank(),
                            panel.background = element_rect(fill = "white", colour = "black"),
                            panel.spacing=unit(0,"cm"))+
        coord_fixed()+labs(fill="log10 FC")
      ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binless_differences.png"),width=15,height=5)
    }
  }
  
  #get statistics
  sub="TBX3_1.5M"
  load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base10k_dfuse",5,"_qmin_",0.01,".RData"))
  resolution=10000
  compnames=data.table(t(matrix(cs@experiments[,name][c(c(1,2),c(1,3),c(2,4),c(3,4))],nrow=2)))
  
  #jacard for binned interactions
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSsig")
  #mat=mat[unclass(bin2)-unclass(bin1)>10]
  #number of significant interactions, and percentage
  mat[,.(num=sum(is.significant),pc=sum(is.significant)/.N),by=name]
  #number of significant interactions
  compnames=data.table(t(matrix(cs@experiments[,name][c(c(1,2),c(1,3),c(2,4),c(3,4))],nrow=2)))
  foreach (nref=compnames[,V1],n=compnames[,V2],.combine=rbind) %do% {
    refmat=mat[name==nref,.(bin1,bin2,ref.is=is.significant)]
    merge(mat[name==n,.(name,bin1,bin2,is=is.significant)],refmat,by=c("bin1","bin2"))[
      ,.(ref=nref,N=sum(is),refN=sum(ref.is),inter=sum(is==T&ref.is==T),union=sum(is==T|ref.is==T),
         jacard=sum(is==T&ref.is==T)/sum(is==T|ref.is==T)),by=name]
  }
  
  #jacard for binless interactions
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSbsig")
  mat[,valid:=.SD[begin2-begin1>10000,.N>4],by=c("patchno","name")]
  mat[,is.maximum:=ifelse(valid==T,is.maximum,F)]
  #number of significant interactions, and percentage
  mat[,.(num=sum(is.maximum),pc=sum(is.maximum)/.N),by=name]
  #number of significant interactions
  compnames=data.table(t(matrix(cs@experiments[,name][c(c(1,2),c(1,3),c(2,4),c(3,4))],nrow=2)))
  foreach (nref=compnames[,V1],n=compnames[,V2],.combine=rbind) %do% {
    refmat=mat[name==nref,.(bin1,bin2,ref.is=is.maximum)]
    merge(mat[name==n,.(name,bin1,bin2,is=is.maximum)],refmat,by=c("bin1","bin2"))[
      ,.(ref=nref,N=sum(is),refN=sum(ref.is),inter=sum(is==T&ref.is==T),union=sum(is==T|ref.is==T),
         jacard=sum(is==T&ref.is==T)/sum(is==T|ref.is==T)),by=name]
  }
  
  #number of overlaps for binless interactions
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSbsig")
  mat[,valid:=.SD[begin2-begin1>10000,.N>4],by=c("patchno","name")]
  mat[,is.maximum:=ifelse(valid==T,is.maximum,F)]
  #number of significant interactions, and percentage
  mat[,.(ismax=is.maximum[1]),by=c("patchno","name")][,.(nmax=sum(ismax),ntot=.N,max.pc=sum(ismax)/.N*100),by=name]
  #number of overlaps with reference calls
  compnames=data.table(t(matrix(cs@experiments[,name][c(c(1,2),c(1,3),c(2,4),c(3,4))],nrow=2)))
  foreach (nref=compnames[,V1],n=compnames[,V2],.combine=rbind) %do% {
    refmat=mat[name==nref,.(bin1,bin2,ref.is=is.maximum)]
    mg=merge(mat[name==n,.(name,bin1,bin2,patchno,is=is.maximum)],refmat,by=c("bin1","bin2"))
    mg[is==T][,.(ov=any(ref.is)),by=c("name","patchno")][,.(ref=nref,n.ov=sum(ov),ntot=.N,ov.pc=100*sum(ov)/.N),by=name]
  }
  
  
}
