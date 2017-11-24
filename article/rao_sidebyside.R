library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

args=commandArgs(trailingOnly=TRUE)
sub=args[1]
ncores=1
base.res=5000

setwd("/home/yannick/simulations/cs_norm")

fit=T
restart=F
bin=T
plot=T

if (fit == T) {
  if (restart == F) {
    load(paste0("data/rao_HiCall_GM12878_",sub,"_csdata.RData"))
    csd1=csd
    load(paste0("data/rao_HiCall_IMR90_",sub,"_csdata.RData"))
    csd2=csd
    cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
    cs <- normalize_binless(cs, restart=F, ncores=ncores, base.res=base.res)
  } else {
    load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
    if (!has_converged(cs)) cs <- normalize_binless(cs, restart=T, ngibbs = 10, ncores=ncores)
  }
  save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
}

if (bin == T) {
  load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  # for (resolution in c(base.res)) {
  #   cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  #   cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  #   cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  #   ref=cs@experiments[1,name]
  #   cs=detect_binned_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
  #   cs=detect_binless_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
  #   save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  # }
  for (resolution in c(10000,20000,50000)) {
    cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
    cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  }
  save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
}

if (plot == T) {
  #load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  for (resolution in c(base.res)) {
    #side-by-side at 5k: observed
    mat=get_binless_interactions(cs)
    plot_binless_matrix(mat, upper="observed", lower="observed")
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_observed.pdf"),width=10,height=5)
    
    #side-by-side at 5k: normalized
    mat=get_binned_interactions(cs, resolution=resolution, group="all")
    mat[,is.significant:=prob.gt.expected>=1-1e-12]
    plot_binless_matrix(mat, upper="normalized", lower="normalized") + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_binned_normalized.pdf"),width=10,height=5)
    
    #side-by-side at 5k: binned differences
    mat=get_binned_differences(cs, resolution=resolution, group="all", ref=cs@experiments[1,name])
    mat[,direction:=ordered(direction,c("enriched","depleted"))]
    plot_binless_difference_matrix(mat) + geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+facet_wrap(~name,nrow = 1)
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_binned_difference.pdf"),width=10,height=5)
    
    #side-by-side at 5k: binless signal
    mat=get_binless_interactions(cs)
    plot_binless_signal_matrix(mat)
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_binless_signal.pdf"),width=10,height=5)
    
    #side-by-side at 5k: binless
    mat=get_binless_interactions(cs)
    plot_binless_matrix(mat, upper="binless", lower="binless")
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_binless.pdf"),width=10,height=5)
    
    #side-by-side at 5k: binless and observed
    mat=get_binless_interactions(cs)
    plot_binless_matrix(mat, upper="binless", lower="observed")
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_binless_observed.pdf"),width=10,height=5)
    
    #side-by-side at 5k: binless difference
    mat=get_binless_differences(cs,ref=cs@experiments[1,name])
    plot_binless_difference_matrix(mat)
    ggsave(filename=paste0("images/rao_",sub,"_GMvsIMR90_base5at5_binless_differences.pdf"),width=8,height=5)
  }
}

if (F) {
  
  subs=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "SEMA3Cext_1.5M", "FOXP1_1.3M", "FOXP1ext_2.3M",
         "TBX3_1.5M", "TBX3ext_2M", "Comparison_1.7M", "Talk_2M", "ADAMTS1_2.3M", "ADAMTS1ext_2.8M")
  registerDoParallel(cores=10)
  info = foreach(sub=subs,.combine=rbind,.errorhandling="remove") %dopar% {
    load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
    data.table(ori=sub,has.converged=binless:::has_converged(cs),#run=run,
               runtime=cs@diagnostics$params[,sum(runtime)],nsteps=cs@diagnostics$params[,max(step)],
               name=cs@experiments[,name],
               lambda1=cs@par$lambda1,lambda2=cs@par$lambda2,eCprime=cs@par$eCprime,lambda_diag=cs@par$lambda_diag,
               binned=length(cs@groups)>0)#,
  }
  info
  
  
  #generate bw plots
  #load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  mat = foreach(resolution=c(5000,10000,20000,50000), .combine=rbind, .errorhandling="remove") %do% {
    mat=get_binned_interactions(cs, resolution=resolution, group="all")[name==name[1]]
    mat[,resolution:=resolution]
  }
  melted.mat=melt(mat[,.(resolution,name,begin1,begin2,observed,decay=decaymat,genomic=biasmat,
                         normalized,normalized.sd=background.sd, signal,signal.sd)],
                  id.vars=c("name","begin1","begin2","resolution"))
  melted.mat[,value:=pmin(value,quantile(value,0.98,na.rm=T)),by=c("variable","name","resolution")]
  melted.mat[,value:=(value-min(value,na.rm=T))/(max(value,na.rm=T)-min(value,na.rm=T)),by=c("variable","name","resolution")]
  #
  ggplot(melted.mat)+geom_raster(aes(begin1,begin2,fill=value))+theme_minimal()+
    geom_raster(aes(begin2,begin1,fill=value))+facet_grid(resolution~variable)+scale_fill_gradient(low="white",high="black",na.value="white")+
    coord_fixed()+theme(panel.background=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="genomic coordinate",y="genomic coordinate")+guides(fill=F)
  ggsave(filename=paste0("images/rao_",sub,"_4resolutions_base5_binned.pdf"),width=32,height=20)
}
