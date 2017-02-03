library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)

#args=commandArgs(trailingOnly=TRUE)
sub="SELP_150k"
bpk=3
type="perf"

setwd("/scratch/workspace/csnorm")

#merge datasets
load(paste0("/scratch/workspace/csnorm_data/data/rao_HiCall_GM12878_",sub,"_csdata.RData"))
csd1=csd
load(paste0("/scratch/workspace/csnorm_data/data/rao_HiCall_IMR90_",sub,"_csdata.RData"))
csd2=csd
cs_stan=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
cs_r=merge_cs_norm_datasets(list(csd2), different.decays="none")
cs_stan=merge_cs_norm_datasets(list(csd2), different.decays="none")

#look at these objects
csd1
csd2
cs

#normalize using approximation
cs_stan = run_gauss(cs_stan, bf_per_kb=bpk, bf_per_decade=10, bins_per_bf=10, ngibbs = 10, iter=100000, init_alpha=1e-7,
                 ncounts = 1000000, type=type, fit_model="stan")
cs_r = run_gauss(cs_r, bf_per_kb=bpk, bf_per_decade=10, bins_per_bf=10, ngibbs = 10, iter=100000, init_alpha=1e-7,
               ncounts = 1000000, type=type, fit_model="exact")
save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_gauss_bpk",bpk,".RData"))

#look at the following objects
cs
a=plot_diagnostics(cs)
a[1]
a[2]
cs_r@par$biases
cs@par$decay

#plot biases and decay
ggplot(cs_r@par$biases)+geom_pointrange(aes(pos,etahat,ymin=etahat-std,ymax=etahat+std,colour=cat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_grid(name ~ cat)#+
ggplot(cs@par$decay)+geom_line(aes(distance,kappa))+
  geom_pointrange(aes(distance,kappahat,ymin=kappahat-std,ymax=kappahat+std), alpha=0.1)+
  facet_wrap(~name,scales = "free")+scale_x_log10()
ggplot(data.table(r=cs_r@par$biases[,etahat],stan=cs_stan@par$biases[,etahat],cat=cs_stan@par$biases[,cat]))+
  geom_point(aes(r,stan))+facet_wrap(~cat)+stat_function(fun=identity)

cs_r@par$alpha
cs_stan@par$alpha


#bin at 20kb
resolution=5000
ncores=30
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ice=100, ncores=ncores)
mat=get_matrices(cs, resolution=resolution, group="all")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(signal)))+
  geom_raster(aes(begin2,begin1,fill=log(signal)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)


#detect significant interactions
cs=detect_interactions(cs, resolution=resolution, group="all", threshold=0.95, ncores=ncores)
mat=get_interactions(cs, type="interactions", resolution=resolution, group="all", threshold=0.95, ref="expected")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)



#detect significant differences
cs=detect_differences(cs, resolution=resolution, group="all", threshold=0.95, ncores=ncores, ref="GM MboI 1")
save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_gauss_bpk",bpk,".RData"))
mat=get_interactions(cs, type="differences", resolution=resolution, group="all", threshold=0.95, ref="GM MboI 1")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)



