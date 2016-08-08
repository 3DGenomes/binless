library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

### read caulobacter dataset and generate with different sampling depths

a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                  skip="SRR",nrows=1000000)
csd=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                     "data/caulo_NcoI_all10M", "WT", "1", skip="SRR", circularize=4042929, dangling.L=c(0,3,5), nrows=10000000,
                     dangling.R=c(3,0,-2), maxlen=600, save.data=F)
cs=merge_cs_norm_datasets(list(csd))
save(cs, file="data/caulo_NcoI_all10M_csnorm.RData")



### normalize different datasets
samplings=c("1M","5M","10M")

registerDoParallel(cores=30)
foreach (sampling=c(samplings,"")) %dopar% {
  load(paste0("data/caulo_NcoI_all",sampling,"_csnorm.RData"))
  bf_per_kb=0.25
  cs = run_simplified(cs, design=NULL, bf_per_kb=bf_per_kb, bf_per_decade=5, bins_per_bf=10, groups=10,
                      lambdas=c(0.1,1,10), ngibbs = 1, iter=100000, ncores=10)
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_all",sampling,"_csnorm_optimized.RData"))
}


#outputs
outputs = foreach (j=c(samplings,""),.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_all",j,"_csnorm_optimized.RData"))
  data.table(dset=j,
             out=tail(cs@par$output,1), runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$value)
}
setkey(outputs,logp)
outputs


### generate plots
prefix="NcoI_all"
samplings=c("1M","5M","10M")

#lFC
lFC = foreach (i=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=i,lFC)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
lFC=rbind(lFC,get_cs_binned(cs,1,"CS")[,.(dset="all",lFC)])
lFC[,dset:=ordered(dset,levels=c("all",samplings))]
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=i,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
mat=rbind(mat,get_cs_binned(cs,1,"CS")[,.(dset="all",begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)])
mat[,dset:=ordered(dset,levels=c("all",samplings))]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black", na.value="white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_normalized.png"), width=10, height=10)

#nu and delta correlation
nu = foreach (i=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  data.table(pos=cs@biases[,pos],log_nu=cs@par$log_nu,log_delta=cs@par$log_delta,dset=i,key="pos")
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
nuref=data.table(pos=cs@biases[,pos],log_nu_ref=cs@par$log_nu,log_delta_ref=cs@par$log_delta,key="pos")
nu=nuref[nu]
nu[,dset:=ordered(dset,levels=samplings)]
ggplot(rbind(nu[,.(var="log_nu",correlation=cor(log_nu_ref,log_nu)),by=dset],
             nu[,.(var="log_delta",correlation=cor(log_delta_ref,log_delta)),by=dset],
             data.table(dset=c("all","all"),var=c("log_nu","log_delta"),correlation=c(1,1))),
       aes(dset,correlation,colour=var))+geom_point()+ylim(0,1)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu_delta_correlation.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=samplings,.combine=rbind,.export="data.table") %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  cs@par$decay[,.(dist,decay,dset=i)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
decay=rbind(decay,cs@par$decay[,.(dist,decay,dset="all")])
decay[,dset:=ordered(dset,levels=c("all",samplings))]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_decay.png"), width=10, height=7.5)







