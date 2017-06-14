library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(rstan)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

sub="SEMA3C_1M"

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


load("data/rao_HiCall_SEMA3C_1M_csnorm_optimized_base10k_bpk30_dfuse20_cpp_poster.RData")
sub="SEMA3C_1M"

#various matrices for GM MboI 1
for (resolution in c(5000,20000)) {
  #observed
  mat=get_matrices(cs, resolution=resolution, group="all")[name=="GM MboI 1"&bin2!=bin2[.N]]
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=log(observed)))+
    geom_raster(aes(begin2,begin1,fill=log(observed)))+coord_fixed()+
    scale_fill_gradient(low="white", high="black",na.value="white")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
    theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/rao_HiCall_GM12878_",sub,"_observed_",resolution/1000,"kb.pdf"), width=10,height=9)
  #
  mat=get_interactions(cs, resolution=resolution, group="all", type="CSbsig")[name=="GM MboI 1"&bin2!=bin2[.N]]
  #phihat
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=(phihat)))+
    geom_raster(aes(begin2,begin1,fill=(phihat)))+coord_fixed()+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"))+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
    theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/rao_HiCall_GM12878_",sub,"_phihat_",resolution/1000,"kb.pdf"), width=10,height=9)
  #phihat.sd
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=sqrt(1/weight)))+
    geom_raster(aes(begin2,begin1,fill=sqrt(1/weight)))+coord_fixed()+
    scale_fill_gradient(low="white",high="black", na.value="black")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
    theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/rao_HiCall_GM12878_",sub,"_phihat_sd_",resolution/1000,"kb.pdf"), width=10,height=9)
  #binless
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=(phi)))+
    geom_raster(aes(begin2,begin1,fill=(phi)))+coord_fixed()+
    scale_fill_gradient2(low=muted("blue"),high=muted("red"))+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
    theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/rao_HiCall_GM12878_",sub,"_binless_",resolution/1000,"kb.pdf"), width=10,height=9)
}

#Interactions and observed matrix for IMR90 at 5kb
resolution=5000
#observed
mat=get_matrices(cs, resolution=resolution, group="all")[name!="GM MboI 1"&bin2!=bin2[.N]]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(observed)))+coord_fixed()+
  scale_fill_gradient(low="white", high="black",na.value="white")+
  scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
  theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.spacing=unit(0,"cm"))
ggsave(filename=paste0("images/rao_HiCall_IMR90_",sub,"_observed_",resolution/1000,"kb.pdf"), width=10,height=9)
#
mat=get_interactions(cs, resolution=resolution, group="all", type="CSbsig")[name!="GM MboI 1"&bin2!=bin2[.N]]
#binless
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=(phi)))+
  geom_raster(aes(begin2,begin1,fill=(phi)))+coord_fixed()+
  scale_fill_gradient2(low=muted("blue"),high=muted("red"))+
  scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
  theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.spacing=unit(0,"cm"))
ggsave(filename=paste0("images/rao_HiCall_IMR90_",sub,"_binless_",resolution/1000,"kb.pdf"), width=10,height=9)

#Binless differences
mat=get_interactions(cs, resolution=resolution, group="all", type="CSbdiff", ref=cs@experiments[1,name])[bin2!=bin2[.N]]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=(delta)))+
  geom_raster(aes(begin2,begin1,fill=(delta)))+coord_fixed()+
  scale_fill_gradient2(low=muted("blue"),high=muted("red"))+
  scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
  theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.spacing=unit(0,"cm"))
ggsave(filename=paste0("images/rao_HiCall_IMR90vsGM_",sub,"_diff_",resolution/1000,"kb.pdf"), width=10,height=9)
#diff reference
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=pmin(5,pmax(-5,phi.ref))))+
  geom_raster(aes(begin2,begin1,fill=pmin(5,pmax(-5,phi.ref))))+coord_fixed()+
  scale_fill_gradient2(low=muted("blue"),high=muted("red"))+
  scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
  theme_void()+ theme(legend.title=element_blank(), axis.title=element_blank(),
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.spacing=unit(0,"cm"))
ggsave(filename=paste0("images/rao_HiCall_IMR90vsGM_",sub,"_ref_",resolution/1000,"kb.pdf"), width=10,height=9)
