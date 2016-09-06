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

#normalize with gibbs sampler
registerDoParallel(cores=2)
foreach (replicate=c("BglIIr1","BglIIr2")) %dopar% {
  load(paste0("data/caulo_",replicate,"_all_csnorm.RData"))
  cs = run_simplified(cs, design=NULL, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=c(0.1,1,10),
                   ngibbs = 1, iter=10000, ncores=30)
  save(cs, file=paste0("data/caulo_",replicate,"_all_csnorm_optimized.RData"))
  cs@binned=list()
  cs=postprocess(cs, resolution=10000, ncores=30, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_",replicate,"_all_csnorm_optimized.RData"))
}



### generate plots
prefix="BglII_all"
fnames=c("data/caulo_NcoI_all_csnorm_optimized.RData",
         "data/caulo_BglIIr1_all_csnorm_optimized.RData",
         "data/caulo_BglIIr2_all_csnorm_optimized.RData")
dsets=c("NcoI", "BglIIr1","BglIIr2")

#outputs
outputs = foreach (i=dsets,j=fnames,.combine=rbind) %do% {
  load(j)
  data.table(dset=i,
             out=tail(cs@par$output,1), runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$value)
}
outputs


#lFC
lFC = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=j,lFC,is.interaction)]
}
ggplot(lFC)+geom_density(aes(lFC,colour=dset))+xlim(-6,4)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_lFC.png"), width=10, height=7.5)
ggplot(lFC[dset=="BglIIr2"])+geom_density(aes(lFC,y=..count..,fill=is.interaction),position="stack")

#normalized matrices
mat = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=j,begin1,begin2,normalized,lFC,observed, expected, is.interaction,prob.observed.gt.expected)]
}
#normalized
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_normalized.png"), width=15, height=5)
#lFC
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=lFC))+
  geom_raster(aes(begin2,begin1,fill=lFC))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~dset)
#observed / expected
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~dset)

#nu and delta correlation
nu = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,dset=j)]
}
nu[,pbin:=cut(pos,3)]
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=dset))+facet_wrap(~pbin,scales = "free_x", nrow=3)+
  scale_y_continuous(limits = c(0,2))+
  ylab("nu")#+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_nu.png"), width=10, height=7.5)
ggplot(nu[pos>5e5&pos<6e5&dset!="NcoI"])+geom_line(aes(pos,exp(log_nu),colour=dset))+scale_y_continuous(limits = c(0,5))

#diagonal decay
decay = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  cs@par$decay[,.(dist,decay,dset=j)]
}
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_decay.png"), width=10, height=7.5)

#general statistics
setkey(matref,begin1,begin2,end1,end2,bin1,bin2)
foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  mat=get_cs_binned(cs,1,"CS")
  setkey(mat,begin1,begin2,end1,end2,bin1,bin2)
  mat=merge(mat,matref,suffixes=c("",".ref"),all=T)
  mat[is.na(is.interaction),is.interaction:=F]
  mat[is.na(is.interaction.ref),is.interaction.ref:=F]
  mat[,.(dset=j,.N),by=c("is.interaction","is.interaction.ref")]
}


