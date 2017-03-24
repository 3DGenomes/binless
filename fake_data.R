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

counts[distance>=1096.47819614&distance<1122.0184543]


