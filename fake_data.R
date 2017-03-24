library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")

csd=generate_fake_dataset()
save(csd,file="data/fake_replicate1_csdata.RData")
csd=generate_fake_dataset(eC=-5,replicate="2")
save(csd,file="data/fake_replicate2_csdata.RData")
csd=generate_fake_dataset(signal=T,eC=-4.3,condition="KO")
save(csd,file="data/fake_signal_replicate1_csdata.RData")
csd=generate_fake_dataset(signal=T,eC=-6.1,condition="KO",replicate="2")
save(csd,file="data/fake_signal_replicate2_csdata.RData")


load("data/fake_replicate1_csdata.RData")
csd1=csd
load("data/fake_replicate2_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
cs = run_gauss(cs, restart=F, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10, ngibbs = 10,
               iter=100000, init_alpha=1e-7, ncounts = 1000000, type="perf", ncores=30)
save(cs,file="data/fake_csnorm_optimized.RData")

ggplot(cs@par$decay)+geom_point(aes(distance,kappahat),alpha=0.1)+
  geom_line(aes(distance,kappa))+facet_wrap(~ name)+scale_x_log10()+
  geom_line(aes(distance,base_count),colour="red",data=cs@counts[sample(.N,min(.N,100000))])


