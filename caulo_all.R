library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


### read caulobacter dataset and generate with different sampling depths

a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                  skip=0L,nrows=1000000)
a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                  skip=0L,nrows=1000000)
csd1=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                      "data/caulo_NcoI_all", "WT", "1", enzyme="NcoI", circularize=4042929, dangling.L=c(0,3,5),
                      dangling.R=c(3,0,-2), maxlen=600, save.data=T)
csd2=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                      "data/caulo_BglIIr1_all", "WT", "1", enzyme="BglII", circularize=4042929, dangling.L=c(0,3,5),
                      dangling.R=c(3,0,-2), maxlen=600, save.data=T)
csd3=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate2_reads_int.tsv",
                      "data/caulo_BglIIr2_all", "WT", "2", enzyme="BglII", circularize=4042929, dangling.L=c(0,3,5),
                      dangling.R=c(3,0,-2), maxlen=600, save.data=T)


load("data/caulo_NcoI_all_csdata.RData")
csd1=csd
load("data/caulo_BglIIr1_all_csdata.RData")
csd2=csd
load("data/caulo_BglIIr2_all_csdata.RData")
csd3=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3), different.decays="none")
save(cs, file="data/caulo_csnorm.RData")



#normalize with gibbs sampler
load("data/caulo_csnorm.RData")
cs = run_simplified(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=c(0.01,0.1,1,10,100),
               ngibbs = 1, iter=10000, ncores=30)
save(cs, file="data/caulo_csnorm_optimized.RData")
load("data/caulo_csnorm_optimized.RData")
cs=bin_all_datasets(cs, resolution=10000, ncores=30, verbose=T, ice=1, dispersion.type=2)
cs@binned[[2]]@individual=detect_interactions(cs, resolution=10000, type="all", dispersion.type=2, dispersion.fun=sum,
                        threshold=0.95, ncores=30, normal.approx=100)
cs=group_datasets(cs, resolution=10000, dispersion.type=3, type="enzyme", dispersion.fun=sum, ice=1, verbose=T)
mat=detect_interactions(cs, resolution=10000, type="enzyme", dispersion.type=3, dispersion.fun=sum,
                        threshold=0.95, ncores=1, normal.approx=100)
#mat2=detect_differences(mat, ref="NcoI", threshold=0.95, ncores=1, normal.approx=100)
save(cs, file="data/caulo_csnorm_optimized.RData")





### generate plots
prefix="caulo"
fnames=c("data/caulo_wt_csnorm_optimized.RData","data/caulo_ko_csnorm_optimized.RData")
#fnames=c("data/caulo_NcoI_all_csnorm_optimized.RData")
#fnames=c("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda0.01.RData")
dsets=c("caulo wt","caulo ko")


#lFC
lFC = foreach (j=c(1,2,3),.combine=rbind) %do% {
  get_matrices(cs, resolution=10000, type="all", dispersion.type=j, dispersion.fun=NA)[,.(name,disp.type=j,lFC=log2(normalized))]
}
ggplot(lFC)+geom_density(aes(lFC,colour=name))+facet_grid(disp.type~.)
#ggsave(filename=paste0("images/",prefix,"_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (j=c(1,2,3),.combine=rbind) %do% {
  get_matrices(cs, resolution=10000, type="all", dispersion.type=j, dispersion.fun=NA)[
    ,.(name,disp.type=j,begin1,begin2,normalized,normalized.sd,icelike,icelike.sd,observed, expected,is.significant,prob.gt.expected)]
}
#icelike
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=prob.gt.expected<0.5),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_grid(disp.type~name)
#ggsave(filename=paste0("images/",prefix,"_normalized.png"), width=10, height=5)
#normalized with error bars
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized.sd)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_grid(disp.type~name)

#nu and delta correlation
nu = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,dset=j)]
}
nu[,pbin:=cut(pos,3)]
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=dset))+facet_wrap(~pbin,scales = "free_x", nrow=3)+
  #scale_y_continuous(limits = c(0,2))+
  ylab("nu")+scale_y_log10()
#ggsave(filename=paste0("images/",prefix,"_nu.png"), width=10, height=7.5)

#sorted with data points
nu.obs = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  a=data.table(id=cs@counts[,id1], pos=cs@counts[,pos1],
             obs=cs@counts[,contact.down]*exp(-cs@pred$log_mean_cdown),
             dset=j)[sample(.N,min(.N,100000))]
  a=merge(a,cs@biases[,.(id,nu=exp(cs@par$log_nu))],by="id")
  a[,obs:=obs*nu]
  setkey(a,nu,pos)
  a[,rank:=.I]
  a
}
ggplot(nu.obs)+geom_line(aes(rank,nu,colour=dset))+geom_point(aes(rank,obs,colour=dset),alpha=0.01)+
  #scale_y_continuous(limits = c(0,2))+
  ylab("nu")+scale_y_log10()


#diagonal decay
decay.obs = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  data.table(dist=cs@counts[,distance],
             obs=cs@counts[,contact.down]*exp(-cs@pred$log_mean_cdown+cs@pred$log_decay_down),
             decay=exp(cs@pred$log_decay_down),
             dset=j)[sample(.N,min(.N,100000))]
}
decay.obs=decay.obs[dist>1000]
setkey(decay.obs, dist)
ggplot(decay.obs)+geom_line(aes(dist,decay,colour=dset))+geom_point(aes(dist,obs,colour=dset),alpha=0.01)+
  scale_x_log10()+scale_y_log10()
#ggsave(filename=paste0("images/",prefix,"_decay.png"), width=10, height=7.5)

