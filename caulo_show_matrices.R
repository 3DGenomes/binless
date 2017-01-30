library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(rstan)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

#normalize
load(paste0("data/caulo_BglII_all_csdata.RData"))
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=10, ngibbs = 20, iter=100000, init_alpha=1e-7, ncounts = 1000000)
save(cs,file=paste0("data/caulo_BglII_all_csnorm_optimized_gauss.RData"))

load("data/caulo_BglII_all_csnorm_optimized_gauss_bpk1_nooutliers_nofill_perf.RData")

ggplot(cs@diagnostics$params[step>4,.(step,leg,value,out.last)])+
  geom_line(aes(step,value))+geom_point(aes(step,value,colour=out.last))+facet_wrap(~leg, scales = "free")+
  theme(legend.position="bottom")
vars=foreach(var=c("eC","eRJ","eDE","alpha","lambda_iota","lambda_rho"),
             trans=(c("exp","exp","exp",NA,"log","log")),.combine=rbind) %do% get_all_values(cs,var,trans)
ggplot(vars[step>4])+geom_line(aes(step,value))+
  geom_point(aes(step,value,colour=leg))+facet_wrap(~variable, scales = "free_y")



ggplot(cs@par$biases)+geom_pointrange(aes(pos,etahat,ymin=etahat-std,ymax=etahat+std,colour=cat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_grid(name ~ cat)#+xlim(500000,600000)#+ylim(-10,10)

ggplot(cs@par$decay)+geom_pointrange(aes(distance,kappahat,ymin=kappahat-std,ymax=kappahat+std,colour=name),alpha=0.1)+
  geom_line(aes(distance,kappa,group=name))+scale_x_log10()

bts=melt(cs@biases,id.vars=c("name","id","pos"))
cts=rbind(cs@counts[,.(name, id=id1, pos=pos1, contact.R=contact.close, contact.L=contact.far)],
          cs@counts[,.(name, id=id1, pos=pos1, contact.R=contact.down,  contact.L=contact.up)],
          cs@counts[,.(name, id=id2, pos=pos2, contact.R=contact.far,   contact.L=contact.close)],
          cs@counts[,.(name, id=id2, pos=pos2, contact.R=contact.down,  contact.L=contact.up)])
cts=rbind(bts,melt(cts[,.(contact.R=sum(contact.R),contact.L=sum(contact.L)),by=c("name","id","pos")],id.vars=c("name","id","pos")))

ggplot(cts)+geom_histogram(aes(value),bins=100)+facet_grid(name~variable,scales="free")
ggplot(cts[,.(value=sum(value)),by=c("name","id")])+geom_histogram(aes(value),bins=100)+facet_grid(~name,scales="free")
cts=cts[name=="WT BglII 2"]
cts[,.(value=sum(value)),by=c("id","name")][value<100,.N]/cts[,.(value=sum(value)),by=c("id","name")][,.N]

cts[,.(value=sum(value)),by=c("id","name")][,quantile(value,c(0.01,0.99))]

load("data/caulo_BglIIr2_all_csdata_with_data.RData")

csd@data[]

#bin
resolution=20000
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ice=1, ncores=ncores)
cs=detect_interactions(cs, resolution=resolution, group="condition", ncores=ncores)
cs=detect_interactions(cs, resolution=resolution, group="condition", ncores=ncores, binless=T)
cs=detect_differences(cs, resolution=resolution, group="condition", ncores=ncores, ref="WT")
cs=detect_differences(cs, resolution=resolution, group="condition", ncores=ncores, ref="WT", binless=T)

mat=get_matrices(cs, resolution=resolution, group="all")
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_normalized_",resolution/1000,"kb.pdf"), width=10,height=9)
#signal
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())+
  #scale_fill_gradient(low="black", high="white")
  scale_fill_gradient2()
ggsave(filename=paste0("images/caulo_BglII_all_signal_",resolution/1000,"kb.pdf"), width=10,height=9)
#observed
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(observed)))+
  scale_fill_gradient(low="white", high="black",na.value="white")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_observed_",resolution/1000,"kb.pdf"), width=10,height=9)
#expected
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=log(expected)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_expected_",resolution/1000,"kb.pdf"), width=10,height=9)


#Interactions at 20kb
mat=get_interactions(cs, resolution=resolution, group="condition", type="interactions", ref="expected", threshold=0.95)
#observed
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_interactions_observed_",resolution/1000,"kb.pdf"), width=19,height=9)
#interactions  
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_interactions_signal_",resolution/1000,"kb.pdf"), width=19,height=9)


#Differences at 20kb
mat=get_interactions(cs, resolution=resolution, group="condition", type="differences", ref="WT", threshold=0.95)
mat=mat[,.(name,signal,ref.signal,begin1,begin2,direction,is.significant,difference)]
mat=rbind(mat,mat[name==name[1],.(name="WT (reference)",signal=ref.signal,ref.signal,begin1,begin2,direction=NA,is.significant=F,difference=1)])

ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(signal)))+
  geom_raster(aes(begin2,begin1,fill=-log(ref.signal)))+
  #geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_differences_datasets_",resolution/1000,"kb.pdf"), width=20,height=19)

#difference matrix
ggplot(mat)+facet_wrap(~name)+
  geom_raster(aes(begin1,begin2,fill=-log(difference)))+
  geom_raster(aes(begin2,begin1,fill=-log(difference)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none", axis.title=element_blank())
ggsave(filename=paste0("images/caulo_BglII_all_both_diffsig_",resolution/1000,"kb.pdf"), width=10,height=9)

