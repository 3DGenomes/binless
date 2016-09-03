library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)

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


#zoom on a portion of the dataset
load("data/caulo_BglIIr2_all_csdata_with_data.RData")
data=csd@data[re.closest1>=2000000&re.closest1<=2150000&re.closest2>=2000000&re.closest2<=2150000]
cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
csd = new("CSdata", info=csd@info, settings=list(circularize=-1),
          data=data, biases=cs_data$biases, counts=cs_data$counts)
save(csd, file="data/caulo_BglIIr2_150k_csdata_with_data.RData")
csd@data=data.table()
save(csd, file="data/caulo_BglIIr2_150k_csdata.RData")
csd2=csd


load("data/caulo_NcoI_150k_csdata.RData")
csd1=csd
load("data/caulo_BglIIr1_150k_csdata.RData")
csd2=csd
load("data/caulo_BglIIr2_150k_csdata.RData")
csd3=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3), different.decays="all")
cs=merge_cs_norm_datasets(list(csd1))
save(cs, file="data/caulo_150k_csnorm.RData")



#normalize with serial sampler
#registerDoParallel(cores=30)
#foreach (bpk=bf_per_kb) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_150k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  bf_per_decade=5
  dmin=1-0.01
  dmax=150000+0.01
  bpk=0.25
  lambda=1
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases, design=cs@design, lambda=lambda,
    bf_per_kb=bpk, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = cs@counts, design = cs@design, dmin=dmin, dmax=dmax,
    bf_per_kb=bpk, bf_per_decade=bf_per_decade, iter=100000, verbose = F, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  op$par$value=op$value
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bpk, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs=bin_all_datasets(cs, resolution=10000, ncores=10, verbose=T, ice=1)
  cs@binned[[1]]=group_datasets(cs@experiments, cs@binned[[1]], type="enzyme", ice=1, verbose=T)
save(cs, file=paste0("data/caulo_NcoI_150k_bfpkb",bpk,"_lambda",lambda,"_csnorm_optimized_new.RData"))
#}

ggplot(cs@binned[[1]]@mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~groupname)

  load("data/caulo_NcoI_150k_bfpkb0.25_lambda1_csnorm_optimized.RData")
  cst=cs
  load("data/caulo_NcoI_150k_bfpkb0.25_lambda1_csnorm_optimized_new.RData")
  c(cst@par$eRJ,cs@par$eRJ[1])
  c(cst@par$eDE,cs@par$eDE[1])
  c(cst@par$eC,cs@par$eC[1])
  #
  ggplot(data.table(id=1:length(cst@par$beta_nu), old=cst@par$beta_nu, new=cs@par$beta_nu[2,]))+
    geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
  mean(cst@par$init$beta_nu)
  mean(init.op$beta_nu[2,])
  #
  ggplot(data.table(id=1:80, old=cst@par$log_nu, new=cs@par$log_nu[81:161]))+
    geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
  mean(cst@par$log_nu)
  mean(cs@par$log_nu)
  #
  ggplot(data.table(id=1:21, old=cs@par$log_nu[81:101], new=cs@par$log_nu[102:122]))+
    geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
  #
  ggplot() +
    geom_line(aes(dist,old),colour="red",data=data.table(dist=cst@counts[,distance], old=cst@par$log_decay, key="dist"))+
    geom_line(aes(dist,new),colour="green",
              data=data.table(dist=cs@counts[,distance], new=cs@par$log_decay,
                              name=cs@counts[,name], key="dist")[name=="WT NcoI 1"])+scale_x_log10()
  
  mean(cst@par$log_nu)
  mean(cs@par$log_nu)
  


### generate plots
prefix="caulo"
fnames=c("data/caulo_wt_csnorm_optimized.RData","data/caulo_ko_csnorm_optimized.RData")
#fnames=c("data/caulo_NcoI_all_csnorm_optimized.RData")
#fnames=c("data/rao_HiCall_chrX_450k_csnorm_optimized_bfpkb1_lambda0.01.RData")
dsets=c("caulo wt","caulo ko")


#lFC
lFC = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=j,lFC)]
}
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
#ggsave(filename=paste0("images/",prefix,"_lFC.png"), width=10, height=7.5)

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
#ggsave(filename=paste0("images/",prefix,"_normalized.png"), width=10, height=5)
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

