library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

args=commandArgs(trailingOnly=TRUE)
base.res=as.integer(args[1])
sub="SEMA3C_1M"
bpk=30
dfuse=5 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
qmin=0.01
ncores=5
setwd("/home/yannick/simulations/cs_norm")

if (F) {
  sub="Fig1C"
  size="1M"
  csd=read_and_prepare(paste0("/scratch/rao/remapping/hg19_chr21/SRR1658572_hg19_chr21_35M-36M_Fig1C.tsv"),
                         paste0("data/rao_HiC003_GM12878_hg19_",sub,"_",size), "GM", run,
                         enzyme="MboI", name=paste(sub,"replicate H"), circularize=-1, dangling.L=c(0),
                         dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)
}

load(paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csdata.RData"))
cs=merge_cs_norm_datasets(list(csd), different.decays="none", dfuse=dfuse, qmin=qmin)
cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 25, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
               ncounts = 1000000, ncores=ncores, base.res=base.res, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
save(cs,file=paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csnorm_optimized_base",base.res/1000,"k_qmin",qmin,".RData"))


for (resolution in c(5000,10000,20000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  save(cs,file=paste0("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csnorm_optimized_base",base.res/1000,"k_qmin",qmin,".RData"))
}  
  
if (F) {
  load("data/rao_HiC003_GM12878_hg19_Fig1C_1M_csnorm_optimized_base10k_qmin0.01_stripped.RData")
  base.res=10000
  
  resolution=5000
  #observed: 06
  mat=get_matrices(cs, resolution=resolution, group="all")
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=observed))+
    geom_raster(aes(begin2,begin1,fill=observed))+coord_fixed()+
    scale_fill_gradient(low="white", high="black",na.value="white",trans="log10")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_observed_",resolution/1000,"kb.pdf"), width=10,height=9)
  #
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSsig")
  #binned interaction: 06
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=log10(signal)))+
    geom_raster(aes(begin2,begin1,fill=log10(signal)))+coord_fixed()+
    geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm")) + labs(fill="log10 FC")
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_binned_signif_",resolution/1000,"kb.pdf"), width=10,height=9)
  #
  idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
  csg=cs@groups[[idx1]]
  idx2=get_cs_interaction_idx(csg, type="CSbsig", raise=T)
  csi=csg@interactions[[idx2]]
  csi@settings$min.l10FC=0.2
  mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
  mat[,value:=phi]
  mat = binless:::detect_binless_patches(mat, csi@settings)
  mat[,value:=NULL]
  mat[,phi.max:=ifelse(is.maximum==T,NA,phi)]
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi/log(10)))+
    geom_raster(aes(begin2,begin1,fill=phi.max/log(10)))+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="black")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))+
    coord_fixed()+labs(fill="log10 FC")
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_binless_signif_",resolution/1000,"kb.pdf"), width=10,height=9)
  
  #binless with decay
  mat=merge(mat,get_matrices(cs, resolution=resolution, group="all")[,.(name,bin1,bin2,decaymat)],by=c("name","bin1","bin2"))
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi/log(10)+log10(decaymat)))+
    geom_raster(aes(begin2,begin1,fill=phi/log(10)+log10(decaymat)))+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="black")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))+
    coord_fixed()+labs(fill="log10 FC")
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_binless_with_decay_",resolution/1000,"kb.pdf"), width=10,height=9)
  
  
  #all bw matrices
  mat=get_matrices(cs, resolution=resolution, group="all")
  melted.mat=melt(mat[,.(name,begin1,begin2,observed=pmin(observed,200),decay=decaymat,genomic=biasmat,
                         normalized=pmin(normalized,10),normalized.sd=pmin(normalized.sd,5))],
                  id.vars=c("name","begin1","begin2"))
  melted.mat[,value:=value/.SD[begin1==begin2,max(value)],by=c("variable","name")]
  melted.mat[variable=="genomic",value:=(value-min(value))/(1-min(value))]
  #one dataset, all bw matrices
  ggplot(melted.mat)+geom_raster(aes(begin1,begin2,fill=value))+
    geom_raster(aes(begin2,begin1,fill=value))+facet_grid(name~variable)+scale_fill_gradient(low="white",high="black",na.value="white")+
    coord_fixed()+theme(panel.background=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_all_binned_",resolution/1000,"kb.pdf"), width=20,height=5)
  
  mat=get_matrices(cs, resolution=resolution, group="all")
  mat.sig=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
  mat.all.sig=merge(mat[,.(name,bin1,bin2,begin1,begin2,observed,decaymat,normalized,signal)],
                    mat.sig[,.(name,bin1,bin2,phihat,weight,beta,phi,patchno)])
  mat.all=melt(mat.all.sig[,.(name,begin1,begin2,log10_normalized=log10(normalized),log10_weight=log10(weight),
                                    log10_signal=log10(signal),
                                    binless=phi/log(10),binless_with_decay=phi/log(10)+log10(decaymat))],
                     id.vars=c("name","begin1","begin2"))
  mat.all[,value:=pmax(-1,pmin(1,value))]
  #one dataset, all colour matrices
  ggplot(mat.all)+geom_raster(aes(begin1,begin2,fill=value))+
    geom_raster(aes(begin2,begin1,fill=value))+facet_grid(name~variable)+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
    coord_fixed()+theme(panel.background=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_all_binless_",resolution/1000,"kb.pdf"), width=20,height=5)
  
  #binless size distribution
  mat=get_matrices(cs, resolution=resolution, group="all")
  mat.sig=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
  mat.all.sig=merge(mat[,.(name,bin1,bin2,begin1,begin2,observed,decaymat,normalized,signal)],
                    mat.sig[,.(name,bin1,bin2,phihat,weight,beta,phi,patchno)])
  patchdistr=mat.all.sig[,.(size=.N,ncounts=sum(observed)),by=c("name","patchno")]
  patchdistr[,surface.pc:=100*size/sum(size)]
  patchdistr[,cat:=ordered(ifelse(size<=4,"small",ifelse(surface.pc>1,"large","medium")),c("small","medium","large"))]
  ggplot(patchdistr)+geom_histogram(aes(surface.pc,fill=cat),bins=100)+scale_x_log10()+
    labs(x="patch surface (% total)",y="count",fill="patch size")
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_patch_surf_",resolution/1000,"kb.pdf"), width=10,height=9)
  
  #reads per patch
  ggplot(patchdistr)+geom_histogram(aes(ncounts,fill=cat),position="stack",bins=50)+scale_x_log10()+
    labs(x="number of reads per patch",y="number of patches", fill="patch size")#+facet_wrap(~cat)
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_patch_nreads_",resolution/1000,"kb.pdf"), width=10,height=9)
  
  #read density
  ggplot(patchdistr)+geom_histogram(aes(ncounts/surface.pc,fill=cat),position="stack",bins=50)+scale_x_log10()+
    labs(x="read density per patch",y="number of patches", fill="patch size")#+facet_wrap(~cat)
  ggsave(filename=paste0("images/rao_HiC003_GM12878_hg19_Fig1C_1M_base",base.res/1000,"k_patch_density_",resolution/1000,"kb.pdf"), width=10,height=9)
  
  #reads per distance
  mat.all.sig[,patchsize:=.N,by=c("name","patchno")]
  mat.all.sig[,distance:=begin2-begin1]
  ggplot(mat.all.sig[,.(var=mean(patchsize)),by=distance])+geom_point(aes(distance,var))
  ggplot(mat.all.sig[,.(var=max(patchsize)),by=distance])+geom_point(aes(distance,var))
  ggplot(mat.all.sig[,.(var=median(patchsize)),by=distance])+geom_point(aes(distance,var))
  ggplot(mat.all.sig[,.(var=min(patchsize)),by=distance])+geom_point(aes(distance,var))
  
  
}
