library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


load("data/rao_HiCall_chrX_450k_csdata_with_data.RData")
begin=73800000
end=73820000
data=cs@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
csd = new("CSdata", info=cs@info, settings=list(circularize=-1),
          data=data, biases=cs_data$biases, counts=cs_data$counts)
save(csd, file="data/rao_HiCall_chrX_20k_csdata_with_data.RData")
csd@data=data.table()
save(csd, file="data/rao_HiCall_chrX_20k_csdata.RData")
csd2=csd



load("data/rao_HiCall_chrX_20k_csdata.RData")
csd1=csd
cs=merge_cs_norm_datasets(list(csd1))


cs=run_exact(cs, bf_per_kb = 5, bf_per_decade = 5, lambdas = 10**seq(from=-1,to=1,length.out=6), ncores = 30, iter = 100000)
save(cs, file="data/rao_HiCall_chrX_20k_csnorm_optimized_exact_bpk5.RData")

run_simplified_gibbs(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=10, groups=10, lambda=50,
                    ngibbs = 3, iter=10000)
cs = run_simplified(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=-1,to=1,length.out=6),
                    ngibbs = 3, iter=10000, ncores=30)
save(cs, file="data/rao_HiCall_chrX_20k_csnorm_optimized_gibbs_grp10_bpk5.RData")

load("data/rao_HiCall_chrX_20k_csnorm_optimized_gibbs_grp10_bpk5.RData")
load("data/rao_HiCall_chrX_20k_csnorm_optimized_exact_bpk5.RData")
plots.e=check_fit(cs)


### Binning at a given resolution
cs=bin_all_datasets(cs, resolution=1000, ncores=30, verbose=T, ice=1)
mat=get_matrices(cs, resolution=1000, group="all")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

### Interaction calling 
cs=detect_interactions(cs, resolution=5000, group="all", detection.type=1, threshold=0.95, ncores=30)
cs=detect_interactions(cs, resolution=5000, group="all", detection.type=2, threshold=0.95, ncores=30)
mat1=get_interactions(cs, type="interactions", resolution=5000, group="all", detection.type=1,
                      threshold=0.95, ref="expected")
mat2=get_interactions(cs, type="interactions", resolution=5000, group="all", detection.type=2,
                     threshold=0.95, ref="expected")
mat=rbindlist(list(det1=mat1,det2=mat2),use.names=T,idcol="method")
#for example, we can plot the ice-like matrices with highlighted interactions in the upper left corner
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=prob.gt.expected<0.5),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_grid(method~name)




### Grouping of datasets

#Datasets can be grouped (merged) to improve the quality of the matrices and the interaction detection power
#for example, we can group datasets by condition and enzyme. For that purpose, pass group the corresponding groupings.
cs=group_datasets(cs, resolution=20000, group=c("condition", "enzyme"), ice=1, verbose=T)

#Again, you can detect interactions in these grouped datasets, and plot the results
cs=detect_interactions(cs, resolution=20000, group=c("condition", "enzyme"), detection.type=1,
                        threshold=0.95, ncores=30)
mat=get_interactions(cs, type="interactions", resolution=20000, group=c("condition", "enzyme"), detection.type=1,
                     threshold=0.95, ref="expected")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=prob.gt.expected<0.5),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)


### Difference calling

#you can call significant differences. It's just like calling interactions, but you need to specify a reference.
#All other matrices will then be compared to that one.
cs=detect_differences(cs, ref="WT BglII", resolution=20000, group=c("condition", "enzyme"), detection.type=1,
                      threshold=0.95, ncores=30)

#the interactions can be retrieved as usual, with a few changes
# specify type="differences" and ref as given in detect_differences
mat=get_interactions(cs, type="differences", resolution=20000, group=c("condition", "enzyme"), detection.type=1,
                     threshold=0.95, ref="WT BglII")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=`prob.gt.WT BglII`<0.5),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

#You can also plot the matrix of the ratio of one experiment to its reference, along with its error bar
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(ratio)))+
  geom_raster(aes(begin2,begin1,fill=log(ratio.sd)))+
  scale_fill_gradient2(na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

save(cs, file="data/caulo_csnorm_optimized.RData")



### END


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

