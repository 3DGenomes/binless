library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(rstan)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

sub="FOXP1"

#normalize
load(paste0("data/rao_HiCall_",sub,"_csdata.RData"))
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=10, ngibbs = 20, iter=100000, init_alpha=1e-7, ncounts = 1000000)
save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_gauss.RData"))

#bin
for (resolution in c(5000,20000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ice=1, ncores=ncores)
  cs=detect_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_differences(cs, resolution=resolution, group="all", ncores=ncores, ref="GM MboI 1")
}


load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_gauss_bpk3_nofill_perf.RData")

#different matrices for GM MboI 1
for (resolution in c(5000,20000)) {
  mat=get_matrices(cs, resolution=resolution, group="condition")#[name=="GM"]
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
    scale_fill_gradient2()+facet_wrap(~name)
  ggsave(filename=paste0("images/rao_HiCall_",sub,"_signal_",resolution/1000,"kb.pdf"), width=10,height=9)
  #observed
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=log(observed)))+
    geom_raster(aes(begin2,begin1,fill=log(observed)))+
    scale_fill_gradient(low="white", high="black",na.value="white")+
    theme_bw()+theme(legend.position = "none", axis.title=element_blank())
  ggsave(filename=paste0("images/rao_HiCall_",sub,"_observed_",resolution/1000,"kb.pdf"), width=10,height=9)
  #expected
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=log(expected)))+
    geom_raster(aes(begin2,begin1,fill=log(expected)))+
    scale_fill_gradient(low="white", high="black")+
    theme_bw()+theme(legend.position = "none", axis.title=element_blank())
  ggsave(filename=paste0("images/rao_HiCall_",sub,"_expected_",resolution/1000,"kb.pdf"), width=10,height=9)
}

#Interactions at 20kb
resolution=20000
mat=get_interactions(cs, resolution=resolution, group="condition", type="interactions", ref="expected", threshold=0.95)
#observed
ggplot(mat)+facet_grid(~name)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_both_observed_",resolution/1000,"kb.pdf"), width=19,height=9)
#interactions  
ggplot(mat)+facet_grid(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=signal.signif>1),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_both_signal_",resolution/1000,"kb.pdf"), width=19,height=9)

#Binless interactions at 20kb
resolution=20000
mat=get_interactions(cs, resolution=resolution, group="condition", type="binteractions", ref="expected", threshold=-1)
plot_binless_matrix(mat)
ggsave(filename=paste0("images/rao_HiCall_",sub,"_both_bsignal_",resolution/1000,"kb.pdf"), width=19,height=9)


#Differences at 20kb
resolution=20000
mat=get_interactions(cs, resolution=resolution, group="condition", type="differences", ref="GM", threshold=0.95)
ggplot(mat)+facet_grid(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(ref.signal)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_both_differences_",resolution/1000,"kb.pdf"), width=10,height=9)

#difference matrix
ggplot(mat)+facet_grid(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(difference)))+
  geom_raster(aes(begin2,begin1,fill=-log(difference)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/rao_HiCall_",sub,"_both_diffsig_",resolution/1000,"kb.pdf"), width=10,height=9)

#Binless differences at 20kb
resolution=20000
mat=get_interactions(cs, resolution=resolution, group="condition", type="bdifferences", ref="GM", threshold=-1)
plot_binless_matrix(mat, minima=T)
ggsave(filename=paste0("images/rao_HiCall_",sub,"_both_bdiffsig_",resolution/1000,"kb.pdf"), width=10,height=9)
