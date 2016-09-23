library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


a=examine_dataset("/scratch/ralph/HiC/3_Mapped/Bcell_Sox2_10Mb_both_filled_map.tsv", skip=0L,nrows=1000000)
csd1=read_and_prepare("/scratch/ralph/HiC/3_Mapped/Bcell_Sox2_10Mb_both_filled_map.tsv",
                      "data/ralph_Bcell_Sox2_10Mb", "Bcell", "1", enzyme="MboI", circularize=-1, dangling.L=c(-1,0,3,4,8),
                      dangling.R=c(4,3,0,-1,-5), maxlen=800, save.data=T)
a=examine_dataset("/scratch/ralph/HiC/3_Mapped/EScell_Sox2_10Mb_both_filled_map.tsv", skip=0L,nrows=1000000)
csd2=read_and_prepare("/scratch/ralph/HiC/3_Mapped/EScell_Sox2_10Mb_both_filled_map.tsv",
                      "data/ralph_EScell_Sox2_10Mb", "EScell", "1", enzyme="MboI", circularize=-1, dangling.L=c(-1,0,3,4,8),
                      dangling.R=c(4,3,0,-1,-5), maxlen=800, save.data=T)


binned=bin_counts(cs@counts,resolution=10000,b1=34500000,e1=35000000,b2=34500000,e2=35000000)
binned=bin_counts(cs@counts,resolution=50000,b1=33500000,e1=36000000,b2=33500000,e2=36000000)
sfac=binned[name==cs@experiments[1,name],sum(N)]/binned[name==cs@experiments[2,name],sum(N)]
binned[name==name[2],N:=N*sfac]
ggplot(binned)+
  geom_raster(aes(begin1,begin2,fill=log(N)))+
  geom_raster(aes(begin2,begin1,fill=log(N)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name,scales="free")

for (cell in c("ES","B")) {
  load(paste0("data/ralph_",cell,"cell_Sox2_10Mb_csdata_with_data.RData"))
  begin=33000000
  end=36000000
  data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
  csd = new("CSdata", info=csd@info, settings=list(circularize=-1),
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  csd@data=data.table()
  save(csd, file=paste0("data/ralph_",cell,"cell_Sox2_3Mb_csdata.RData"))
}


load("data/ralph_EScell_Sox2_0.5Mb_csdata.RData")
csd1=csd
load("data/ralph_Bcell_Sox2_0.5Mb_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
save(cs,file="data/ralph_Sox2_0.5Mb_csnorm.Rdata")




load("data/ralph_Sox2_0.5Mb_csnorm.Rdata")

init.op <- csnorm:::csnorm_simplified_guess(
  biases = cs@biases, counts = cs@counts, design = cs@design, lambda=1, dmin=cs@settings$dmin, dmax=cs@settings$dmax,
  groups = 10, bf_per_kb = 1, bf_per_decade = 5, iter = 10000)


cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=0,to=2,length.out=6),
                    ngibbs = 3, iter=10000, ncores=30)
save(cs, file="data/ralph_Sox2_0.5Mb_csnorm_optimized_bpk1_grp10.Rdata")

load("data/ralph_Sox2_0.5Mb_csnorm_optimized_gibbs_bpk1_grp10.Rdata")

for (bpk in 1:5) for (groups in c(10,50,100)) {
  load(paste0("data/ralph_Sox2_0.5Mb_csnorm_optimized_gibbs_bpk",bpk,"_grp",groups,".Rdata"))
  cat(paste("bpk",bpk,"groups",groups,"dispersion",cs@par$alpha,"\n"))
}


plots=check_fit(cs)


### Binning at a given resolution
cs=bin_all_datasets(cs, resolution=10000, ncores=30, verbose=T, ice=1)
mat=get_matrices(cs, resolution=1000, group="all")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

### Interaction calling 
cs=detect_interactions(cs, resolution=20000, group="all", threshold=0.95, ncores=30)
mat=get_interactions(cs, type="interactions", resolution=20000, group="all", threshold=0.95, ref="expected")
#for example, we can plot the ice-like matrices with highlighted interactions in the upper left corner
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=prob.gt.expected<0.5),data=mat1[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_grid(~name)

cs=detect_differences(cs, resolution=10000, group="all", detection.type=1, threshold=0.95, ncores=30, ref="EScell MboI 1")
cs=detect_differences(cs, resolution=10000, group="all", detection.type=2, threshold=0.95, ncores=30, ref="EScell MboI 1")
cs=detect_differences(cs, resolution=10000, group="all", detection.type=3, threshold=3, ncores=30, ref="EScell MboI 1")
#
mat1=get_interactions(cs, type="differences", resolution=10000, group="all", detection.type=1, threshold=0.95, ref="EScell MboI 1")
mat2=get_interactions(cs, type="differences", resolution=10000, group="all", detection.type=2, threshold=0.95, ref="EScell MboI 1")
mat3=get_interactions(cs, type="differences", resolution=10000, group="all", detection.type=3, threshold=3, ref="EScell MboI 1")
mat1[,c("ratio","ratio.sd"):=list(NULL,NULL)]
mat=rbindlist(list(det1=mat1,det2=mat2,det3=mat3),use.names=T,idcol="method")
#for example, we can plot the ice-like matrices with highlighted interactions in the upper left corner
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat1[is.significant==T])+
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

