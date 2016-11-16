library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(rstan)

setwd("/home/yannick/simulations/cs_norm")

sub="GM12878_FOXP1_newmat_sigmainf"

#normalize
load(paste0("data/rao_HiCall_",sub,"_csdata.RData"))
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=10, ngibbs = 20, iter=100000, init_alpha=1e-7, ncounts = 1000000)
save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_gauss.RData"))

#bin
resolution=5000
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ice=100, ncores=ncores)
cs=group_datasets(cs, resolution=50000, group="condition")
for (i in 1:6) cs@binned[[i]]@grouped[[1]]@interactions=list()
for (resolution in c(5000,10000,15000,16000,18000,20000)) {
  cs=detect_interactions(cs, resolution=resolution, group="all", ncores=ncores)
}

#normalized matrix
mat=get_matrices(cs, resolution=resolution, group="all")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_normalized_",resolution/1000,"kb.pdf"), width=10,height=9)
#signal
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())+
  #scale_fill_gradient(low="black", high="white")
  scale_fill_gradient2()
ggsave(filename=paste0("images/rao_HiCall_",sub,"_signal_",resolution/1000,"kb.pdf"), width=10,height=9)
#observed
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(observed)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_observed_",resolution/1000,"kb.pdf"), width=10,height=9)
#expected
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(expected)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_expected_",resolution/1000,"kb.pdf"), width=10,height=9)


#interactions
mat=foreach (resolution=c(5000,20000),.combine=rbind) %do% {
  mat=get_interactions(cs, resolution=resolution, group="all", ref="expected", type="interactions", threshold=0.95)
  mat[,resolution:=resolution]
  mat
}
ggplot(mat)+#facet_wrap(~ resolution)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal.signif)),data=mat[is.significant==T])+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())+
  scale_fill_gradient(low="black", high="white")
  #scale_fill_gradient2(na.value="black")
ggsave(filename="foxp1_alpha10_beta15.png",width=15,height=10)

ggplot(mat)+geom_point(aes(signal,signal.signif))+stat_function(fun=identity)+scale_x_log10()+scale_y_log10()
ggplot(mat)+geom_point(aes(signal.sd,signal.signif.sd))+stat_function(fun=identity)+scale_x_log10()+scale_y_log10()
