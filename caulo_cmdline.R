library(csnorm)
#setwd("/home/yannick/simulations/cs_norm")
#load("data/caulo_BglIIr1_500k_csdata.RData")
#csd1=csd
#load("data/caulo_BglIIr2_500k_csdata.RData")
#csd2=csd
#load("data/caulo_BglII_rifampicin_500k_csdata.RData")
#csd3=csd
#cs=merge_cs_norm_datasets(list(csd1,csd2,csd3), different.decays="none")
#save(cs,file="step0.RData")
#cs = run_simplified(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=-2,to=1,length.out=6),
#                    ngibbs = 3, iter=10000, ncores=30)
#save(cs,file="step1.RData")
#cs=bin_all_datasets(cs, resolution=20000, ncores=30, verbose=T, ice=1)
#save(cs,file="step2.RData")
load("step2.RData")
cat("DETECT INTERACTIONS")
cs=detect_interactions(cs, resolution=20000, ncores=30, group="all", detection.type=2)
save(cs,file="step3.RData")
