library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

ncores=5
base.res=5000

load(paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csdata.RData"))
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = normalize_binless(cs, restart=F, ncores=ncores, base.res=base.res)
save(cs,file=paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csnorm_optimized_base",base.res/1000,"k.RData"))


resolution=base.res
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
save(cs,file=paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csnorm_optimized_base",base.res/1000,"k.RData"))

if (F) {
  load(paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csnorm_optimized_base",base.res/1000,"k.RData"))
  
  resolution=5000
  #observed
  mat=get_binned_matrices(cs, resolution=resolution, group="all")
  plot_binless_matrix(mat,upper="observed",lower="observed")
  
  #normalized with binned interactions
  mat=get_binned_interactions(cs, resolution=resolution, group="all")
  plot_binless_matrix(mat, upper="normalized", lower="normalized") + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+facet_wrap(~name,nrow = 1)
  
  #binless signal
  mat=get_binless_interactions(cs)
  plot_binless_signal_matrix(mat)+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_binless_signal.pdf"), width=10,height=9)
  
  #binless matrix
  mat=get_binless_interactions(cs)
  plot_binless_matrix(mat, upper="binless", lower="binless")+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_binless.pdf"), width=10,height=9)
  
  #binless size distribution
  mat.all.sig=get_binless_interactions(cs)
  patchdistr=mat.all.sig[,.(size=.N,ncounts=sum(observed)),by=c("name","patchno")]
  patchdistr[,surface.pc:=100*size/sum(size)]
  patchdistr[,cat:=ordered(ifelse(size<=4,"small",ifelse(surface.pc>1,"large","medium")),c("small","medium","large"))]
  ggplot(patchdistr)+geom_histogram(aes(surface.pc,fill=cat),bins=100)+scale_x_log10()+
    labs(x="patch surface (% total)",y="count",fill="patch size")
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_patch_surf.pdf"), width=10,height=9)
  
  #reads per patch
  ggplot(patchdistr)+geom_histogram(aes(ncounts,fill=cat),position="stack",bins=50)+scale_x_log10()+
    labs(x="number of reads per patch",y="number of patches", fill="patch size")#+facet_wrap(~cat)
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_patch_nreads.pdf"), width=10,height=9)
  
  #read density
  ggplot(patchdistr)+geom_histogram(aes(ncounts/surface.pc,fill=cat),position="stack",bins=50)+scale_x_log10()+
    labs(x="read density per patch",y="number of patches", fill="patch size")#+facet_wrap(~cat)
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_patch_density.pdf"), width=10,height=9)
  
  #patches per distance
  mat.all.sig[,patchsize:=.N,by=c("name","patchno")]
  mat.all.sig[,distance:=begin2-begin1]
  ggplot(mat.all.sig[,.(var=mean(patchsize)),by=distance])+geom_point(aes(distance,var))
  ggplot(mat.all.sig[,.(var=max(patchsize)),by=distance])+geom_point(aes(distance,var))
  ggplot(mat.all.sig[,.(var=median(patchsize)),by=distance])+geom_point(aes(distance,var))
  ggplot(mat.all.sig[,.(var=min(patchsize)),by=distance])+geom_point(aes(distance,var))
  
  
}
