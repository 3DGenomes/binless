library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

sub="SELP_150k"

#load(paste0("data/rao_HiCall_GM12878_",sub,"_csdata.RData"))
#csd1=csd
#load(paste0("data/rao_HiCall_IMR90_",sub,"_csdata.RData"))
#csd2=csd
#cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")

load("data/rao_HiCall_GM12878_SELP_150k_csnorm_optimized_gauss.RData")
load("data/rao_HiCall_GM12878_SELP_150k_csnorm_optimized_exact.RData")
cs@counts=fill_zeros(cs@counts, cs@biases, circularize = cs@settings$circularize, dmin = cs@settings$dmin)

cs@par=list()
cs=run_gauss(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=20, ngibbs = 5,
             iter=100000, init_alpha=1e-7, ncounts = 1000000, type="outer")

#pred=csnorm_predict_all(cs,cs@counts)
pred=csnorm_predict_all_parallel(cs, cs@counts, ncores=30)
#decay mean
pred[,mean(log_decay)]
cs@par$decay[,weighted.mean(kappa,ncounts)-cs@par$eC]
#iota mean
pred[id1==1,mean(log_mean_cfar-cs@par$eC-log_decay)]
mean(cs@par$log_iota)
cs@par$biases[cat=="dangling L",mean(eta)-cs@par$eDE]
#rho mean
pred[id1==1,mean(log_mean_cclose-cs@par$eC-log_decay)]
mean(cs@par$log_rho)
cs@par$biases[cat=="dangling R",mean(eta)-cs@par$eDE]

