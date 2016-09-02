library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

### read rao dataset

a=examine_dataset("/scratch/rao/mapped/HICall_both_filled_map_chrX_73780165-74230165.tsv",
                  skip="SRR",nrows=1000000)
csd=read_and_prepare("/scratch/rao/mapped/HICall_both_filled_map_chrX_73780165-74230165.tsv",
                     "data/rao_HiCall_chrX_450k", "WT", "1", skip="SRR", circularize=-1, dangling.L=c(-1,0,3,4,5,8),
                     dangling.R=c(4,3,0,-1,-2,-5), maxlen=1000, save.data=T)
cs=merge_cs_norm_datasets(list(csd))
save(cs, file="data/rao_HiCall_chrX_450k_csnorm.RData")


### normalize different datasets

bf_per_kb=1
bf_per_decade=5
bins_per_bf=10
groups=10
lambdas=c(0.01,0.1,1)
registerDoParallel(cores=3)
foreach (lambda=lambdas) %dopar% {
  load("data/rao_HiCall_chrX_450k_csnorm.RData")
  cs = run_gibbs(cs, design=NULL, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf,
                 groups=groups, lambda=lambda, ngibbs = 1, iter=10000)
  save(cs, file=paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",lambda,".RData"))
}
foreach (lambda=lambdas) %dopar% {
  load(paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",lambda,".RData"))
  cs@binned=list()
  cs=postprocess(cs, resolution=5000, ncores=30, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",lambda,".RData"))
}


prefix="rao_HiCall_chrX_450k"

#outputs and runtime
outputs = foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",j,".RData"))
  data.table(dset=j,
             out=tail(cs@par$output,1), runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$value)
}
setkey(outputs,logp)
outputs

#nu and delta
nu = foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",j,".RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,dset=j)]
}
ggplot(nu)+geom_line(aes(pos,log_nu))+facet_grid(dset~.)
ggplot(nu[pos>=73800000&pos<=73900000])+geom_line(aes(pos,log_nu))+facet_grid(dset~.)
ggplot(nu)+geom_line(aes(pos,log_delta))+facet_grid(dset~.)
ggplot(nu[pos>=73800000&pos<=73900000])+geom_line(aes(pos,log_delta))+facet_grid(dset~.)


### generate plots
#lFC
lFC = foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",j,".RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=j,lFC)]
}
lFC[,dset:=ordered(dset,levels=lambdas)]
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda",j,".RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=j,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
mat[,dset:=ordered(dset,levels=lambdas)]
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


