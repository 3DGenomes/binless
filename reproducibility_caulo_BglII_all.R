library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")

### read caulobacter dataset and generate with different sampling depths

a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                  skip="SRR",nrows=1000000)
csd=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                     "data/caulo_BglIIr1_all", "WT", "1", skip="SRR", circularize=4042929, dangling.L=c(0,3,5),
                     dangling.R=c(3,0,-2), maxlen=600, save.data=T)
cs=merge_cs_norm_datasets(list(csd))
save(cs, file="data/caulo_BglIIr1_all_csnorm.RData")



### normalize different datasets

load("data/caulo_BglIIr1_all_csnorm.RData")

coverage=4
square.size=500000
bf_per_kb=0.25
cs=run_split_parallel(cs, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                      bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=2000, ncores=30,
                      homogenize=F, outprefix="tmp/test")#, ops.count=ops.count, ops.bias=ops.bias)
#cs=run_split_parallel_recovery(cs, "tmp/test", square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
#                      bf_per_decade=5, distance_bins_per_decade=100, verbose = F, iter=10000, ncores=30,
#                      homogenize=F)
#a=csnorm:::diagnose_biases(cs, outprefix = "tmp/test", coverage = coverage, square.size = square.size)
#a=csnorm:::diagnose_counts(cs, outprefix = "tmp/test", coverage.extradiag=1, square.size = square.size)
cs@pred=csnorm_predict_all(cs, ncores=30)
cs=postprocess(cs, resolution=10000, ncores=30, verbose=F)
cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
save(cs, file="data/caulo_BglIIr2_all_csnorm_optimized_bfpkb025.RData")
#save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))



### generate plots
prefix="BglII_all"
fnames=c("data/caulo_BglIIr1_all_csnorm_optimized.RData",
         "data/caulo_BglIIr1_all_csnorm_optimized_redo.RData",
         "data/caulo_BglIIr1_all_csnorm_optimized_bfpkb025.RData",
         "data/caulo_BglIIr2_all_csnorm_optimized.RData")
dnames=c("BglIIr1","BglIIr1 redo", "BglIIr1 0.25", "BglIIr2")
#lFC
lFC = foreach (i=fnames,d=dnames,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=d,lFC)]
}
lFC[,dset:=ordered(dset,levels=dnames)]
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=fnames,d=dnames,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=d,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
mat[,dset:=ordered(dset,levels=dnames)]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_normalized.png"), width=10, height=7.5)

#nu and delta plots
nu = foreach (i=fnames,d=dnames,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],log_nu=cs@par$log_nu,log_delta=cs@par$log_delta,dset=d,key="pos")
}
nu[,dset:=ordered(dset,levels=dnames)]
ggplot(melt(nu,id.var=c("pos","dset")))+geom_line(aes(pos,exp(value),colour=variable))+facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_nu_delta.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=fnames,d=dnames,.combine=rbind) %do% {
  load(i)
  cs@par$decay[,.(dist,decay,dset=d)]
}
decay[,dset:=ordered(dset,levels=dnames)]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_decay.png"), width=10, height=7.5)


