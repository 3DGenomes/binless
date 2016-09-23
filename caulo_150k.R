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
csd4=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_rifampicin_reads_int.tsv",
                      "data/caulo_BglII_rifampicin_all", "rifampicin", "1", enzyme="BglII", circularize=4042929, dangling.L=c(0,3,5),
                      dangling.R=c(3,0,-2), maxlen=600, save.data=T)


#zoom on a portion of the dataset
begin=250000
end=750000
begin=2000000
end=2150000

for (i in c("BglIIr1","BglIIr2","BglII_rifampicin")) {
  load(paste0("data/caulo_",i,"_all_csdata_with_data.RData"))
  data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
  csd = new("CSdata", info=csd@info, settings=list(circularize=-1),
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  save(csd, file=paste0("data/caulo_",i,"_500k_csdata_with_data.RData"))
  csd@data=data.table()
  save(csd, file=paste0("data/caulo_",i,"_500k_csdata.RData"))
  csd2=csd
}


load("data/caulo_BglIIr1_500k_csdata.RData")
csd1=csd
load("data/caulo_BglIIr2_500k_csdata.RData")
csd2=csd
load("data/caulo_BglII_rifampicin_500k_csdata.RData")
csd3=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3), different.decays="none")
save(cs, file="data/caulo_rif_500k_csnorm.RData")

load("data/caulo_BglII_rif_500k_csnorm_optimized.RData")


#normalize with serial sampler
load("data/caulo_rif_500k_csnorm.RData")
cs=run_exact(cs, bf_per_kb = 0.25, bf_per_decade = 5, lambdas = 10**seq(from=-2,to=1,length.out=6), ncores = 30, iter = 100000)
save(cs, file="data/caulo_rif_500k_csnorm_optimized_exact.RData")

#normalize with gibbs sampler
cs = run_simplified(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=-2,to=1,length.out=6),
                    ngibbs = 5, iter=10000, ncores=30)
save(cs, file="data/caulo_BglII_rif_csnorm_optimized_gibbs3_bpk0.25_grp10.RData")

for (d in 2:length(op$par$beta_diag)) {
  if (abs(op$par$beta_diag[d]-op$par$beta_diag[d-1])<.Machine$double.eps) {
    op$par$beta_diag[d:length(op$par$beta_diag)]=op$par$beta_diag[d:length(op$par$beta_diag)]+.Machine$double.eps
  }
}

op.diag <- csnorm:::csnorm_simplified_decay(
  biases = cs@biases, counts = cs@counts, design=cs@design,
  log_nu = cs@par$log_nu, log_delta = cs@par$log_delta,
  dmin = cs@settings$dmin, dmax = cs@settings$dmax, bf_per_decade = cs@settings$bf_per_decade,
  bins_per_bf = cs@settings$bins_per_bf, groups = cs@settings$groups,
  iter=10000, init=cs@par, dispersion=cs@par$alpha, verbose=T, tol_rel_grad=0)
ggplot()+scale_x_log10()+scale_y_log10()+
  geom_line(aes(dist,decay,colour="ori"),data=cs@par$decay)+
  geom_line(aes(dist,decay,colour="new"),data=op.diag.old$par$decay)+
  geom_line(aes(dist,decay,colour="last"),data=op.diag$par$decay)

op.diag.old=op.diag

load("data/caulo_rif_500k_csnorm_optimized_exact.RData")
load("data/caulo_rif_500k_csnorm_optimized_gibbs.RData")

cs=bin_all_datasets(cs, resolution=20000, ncores=30, verbose=T, ice=1)
cs=detect_interactions(cs, resolution=20000, ncores=30, group="all", threshold=0.95)
cs=detect_differences(cs, resolution=20000, ncores=30, group="all", ref="WT BglII 2", threshold=0.95)

mat=get_interactions(cs, type="interactions", resolution=20000, group="all", ref="expected", threshold=0.95)
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient(na.value = "white",low="white",high="black")+theme(legend.position = "none")+
  facet_wrap(~name)
ggsave(filename = "images/caulo_BglII_rif_500k_interactions_20k.png", width=12, height=4)

int=detect_binless_interactions(cs, group="all", ncores=30)
counts=int$counts
predicted=int$predicted
ggplot() + 
  geom_raster(aes(x=pos1,y=pos2,fill=log(signal)),data=predicted)+
  geom_point(aes(pos2,pos1,colour=log(signal)),data=counts)+
  facet_wrap(~groupname)+
  scale_colour_gradient2(limits=c(-1,1))+scale_fill_gradient2(limits=c(-1,1))
ggsave(filename = "images/caulo_BglII_rif_500k_interactions_binless.png", width=12, height=4)

mat1=get_interactions(cs, type="differences", resolution=20000, group="all", ref="WT BglII 2", detection.type=1, threshold=0.95)
mat2=get_interactions(cs, type="differences", resolution=20000, group="all", ref="WT BglII 2", detection.type=2, threshold=0.95)
mat3=get_interactions(cs, type="differences", resolution=20000, group="all", ref="WT BglII 2", detection.type=3, threshold=3)
setnames(mat3,"logK","prob.gt.WT BglII 2")
mat1[,c("ratio","ratio.sd"):=list(NULL,NULL)]
mat=rbindlist(list(detect1=mat1,detect2=mat2,detect3=mat3),use.names=T, idcol="method")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=`prob.gt.WT BglII 2`<0.5),data=mat[is.significant==T])+
  scale_fill_gradient(na.value = "white",low="white",high="black")+theme(legend.position = "none")+
  facet_grid(method~name)

mat2=get_predicted_subset(cs)
setkey(mat,norm,name,distance)
ggplot(mat)+geom_line(aes(distance,exp(log_decay)),colour="red")+geom_point(aes(distance,contact.close/exp(log_mean_cclose-log_decay)),alpha=0.1)+
  scale_y_log10()+scale_x_log10(limits=c(1e3,NA))+facet_grid(norm~name)
ggplot(mat)+geom_line(aes(distance,exp(log_decay),colour=norm))+
  scale_y_log10()+scale_x_log10(limits=c(1e3,NA))+facet_grid(~name)
ggplot(mat[distance<300000])+geom_line(aes(distance,exp(log_decay)),colour="red")+geom_point(aes(distance,contact.close/exp(log_mean_cclose-log_decay)),alpha=0.1)+
  scale_y_log10()+facet_grid(norm~name)




counts=fill_zeros(cs@counts,cs@biases)
setkey(counts,name,id1,id2)
pred=csnorm:::csnorm_predict_all(cs, counts, verbose=F)
bcounts=pred[name=="WT BglII 1"&pos1>=262840&pos2>=402840&pos1<=282840&pos2<=422840] #prob ~ 0.5
bcounts=pred[name=="WT BglII 2"&pos1>=322840&pos2>=362840&pos1<=342840&pos2<=382840] #prob ~ 0
data=list(N=4*bcounts[,.N], observed=bcounts[,c(contact.close,contact.down,contact.far,contact.up)],
          log_expected=bcounts[,c(log_mean_cclose,log_mean_cdown,log_mean_cfar,log_mean_cup)], alpha=cs@par$alpha)
c(mean(data$observed),mean(exp(data$log_expected)))
fit=stan("detection.stan",data=data,init=0)
my_sso <- launch_shinystan(fit)
sm=stan_model("detection.stan")
op=optimizing(sm, data=data, as_vector=F, hessian=T, iter=10000, verbose=F, init=0)
#p(log_s>0)?
pnorm(0,mean=op$par$log_s,sd=1/sqrt(-op$hessian[1,1]),lower.tail=F,log.p=F)


cs=detect_differences(cs, resolution=20000, ncores=30, dispersion.type=3, type="all", dispersion.fun=NA,
                      ref="WT BglII 2", threshold=0.95, normal.approx=100)
cs=group_datasets(cs, resolution=20000, dispersion.type=3, type="enzyme", dispersion.fun=sum, ice=1, verbose=T)
cs=detect_differences(cs, resolution=20000, ncores=30, dispersion.type=3, type="enzyme", dispersion.fun=sum,
                      ref="BglII", threshold=0.95, normal.approx=100)

get_matrices(cs, resolution=20000, dispersion.type=3, type="enzyme", dispersion.fun=sum)
get_interactions(cs, interaction.type="interactions", resolution=20000, type="all", ref="expected", dispersion.type=3, dispersion.fun=NA,
                 threshold=0.95, normal.approx=100)
get_interactions(cs, interaction.type="differences", resolution=20000, type="all", ref="WT BglII 2", dispersion.type=3, dispersion.fun=NA,
                 threshold=0.95, normal.approx=100)
get_interactions(cs, interaction.type="differences", resolution=20000, type="enzyme", ref="BglII", dispersion.type=3, dispersion.fun=sum,
                 threshold=0.95, normal.approx=100)







ggplot(mat[name!="WT BglII 2"])+
  geom_raster(aes(begin1,begin2,fill=log(ratio)))+
  geom_raster(aes(begin2,begin1,fill=log(ratio)))+
  geom_point(aes(begin1,begin2,colour=`prob.gt.WT BglII 2`<0.5),data=mat[is.significant==T])+
  scale_fill_gradient2(na.value = "white")+theme(legend.position = "none")+
  facet_grid(name~disp.type)



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

