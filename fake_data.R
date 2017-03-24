library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")

csd=generate_fake_dataset(signal=F,eC=-4.4,replicate="1",condition="WT")
save(csd,file="data/fake_replicate1_csdata.RData")
biases.ref=csd@biases
csd=generate_fake_dataset(signal=F,biases.ref=biases.ref,eC=-5,replicate="2",condition="WT")
save(csd,file="data/fake_replicate2_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.3,replicate="1",condition="KO")
save(csd,file="data/fake_signal_replicate1_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-6.1,replicate="2",condition="KO")
save(csd,file="data/fake_signal_replicate2_csdata.RData")

csd@counts[,bin1:=round(pos1/10000)]
csd@counts[,bin2:=round(pos2/10000)]
binned=csd@counts[,.(count=sum(contact.far+contact.close+contact.up+contact.down)),by=c("bin1","bin2")]
ggplot(binned)+geom_raster(aes(bin1,bin2,fill=log(count)))


load("data/fake_replicate1_csdata.RData")
csd1=csd
load("data/fake_replicate2_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
cs = run_gauss(cs, restart=F, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10, ngibbs = 20,
               iter=100000, init_alpha=1e-7, ncounts = 1000000, type="perf", ncores=30)
save(cs,file="data/fake_csnorm_optimized.RData")

plot_diagnostics(cs)

ggplot(cs@par$decay)+geom_point(aes(distance,kappahat),alpha=0.1)+
  geom_line(aes(distance,kappa))+facet_wrap(~ name)+scale_x_log10()+
  geom_line(aes(distance,base_count),colour="red",data=cs@counts[sample(.N,min(.N,100000))])

ggplot(cs@par$biases[cat=="dangling L"])+geom_point(aes(pos,etahat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_wrap(~ name)+xlim(0,100000)+
  geom_line(aes(pos,true_log_mean_DL),colour="red",data=cs@biases)

ggplot(cs@par$biases[cat=="dangling R"])+geom_point(aes(pos,etahat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_wrap(~ name)+xlim(0,100000)+
  geom_line(aes(pos,true_log_mean_DR),colour="red",data=cs@biases)

resolution=10000
ncores=30
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
cs=detect_binned_interactions(cs, resolution=resolution, group="all", threshold=0.95, ncores=ncores)
save(cs,file="data/fake_csnorm_optimized.RData")
mat=get_interactions(cs, type="interactions", resolution=resolution, group="all", threshold=0.95, ref="expected")
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=normalized))+facet_wrap(~name)+
  scale_fill_gradient(high="red", low="white", na.value = "white")
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=signal))+facet_wrap(~name)+
  scale_fill_gradient(high="black", low="white", na.value = "white")+
  geom_point(aes(bin1,bin2,colour=direction),data=mat[is.significant==T])
mat[,.N,by=is.significant]


cs=detect_binned_differences(cs, resolution=resolution, group="all", threshold=0.95,
                             ncores=ncores, ref=as.character(cs@experiments[1,name]))
save(cs,file="data/fake_csnorm_optimized.RData")
mat=get_interactions(cs, type="differences", resolution=resolution, group="all", threshold=0.95,
                     ref=as.character(cs@experiments[1,name]))
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=difference))+facet_wrap(~name)+
  scale_fill_gradient(high="black", low="white", na.value = "white")+
  geom_point(aes(bin1,bin2,colour=direction),data=mat[is.significant==T])


cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=30, niter=20, fit.decay=T)
save(cs,file="data/fake_csnorm_optimized.RData")





cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores,
                              ref=as.character(cs@experiments[1,name]), niter=2)
save(cs, file=fname)


