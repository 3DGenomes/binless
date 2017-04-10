library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")


#sub="Peak1_450k"
sub="FOXP1_1.3M"

#plot convergence
a=foreach (BICtype=1:4,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_GM12878_",sub,"_csnorm_optimized_base20k_BIC",BICtype,".RData"))
  cs@diagnostics$params[,.(step,leg,value,BIC=BICtype)]
}
ggplot(a)+geom_point(aes(step,value))+facet_wrap(BIC~leg,scales="free")

#bin data
a=foreach (BICtype=1:4,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_GM12878_",sub,"_csnorm_optimized_base20k_BIC",BICtype,".RData"))
  cs=bin_all_datasets(cs, resolution=cs@settings$base.res, verbose=T, ncores=ncores)
  save(cs,file=paste0("data/rao_HiCall_GM12878_",sub,"_csnorm_optimized_base20k_BIC",BICtype,".RData"))
}

#plot binned matrix
a=foreach (BICtype=1:4,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_GM12878_",sub,"_csnorm_optimized_base20k_BIC",BICtype,".RData"))
  get_matrices(cs, resolution=cs@settings$base.res, group="all")[,.(name,bin1,bin2,observed,expected,BIC=BICtype)]
}
ggplot(a)+geom_raster(aes(bin1,bin2,fill=observed))+geom_raster(aes(bin2,bin1,fill=expected))+
  facet_wrap(~BIC)+scale_fill_gradient(high=muted("red"), low="white", na.value = "white")

#plot last signal
a=foreach (BICtype=1:4,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_GM12878_",sub,"_csnorm_optimized_base20k_BIC",BICtype,".RData"))
  cs@par$signal[,.(name,bin1,bin2,phi,BIC=BICtype)]
}
ggplot(a)+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+
  facet_wrap(~BIC)+scale_fill_gradient(high=muted("red"), low="white", na.value = "white")



