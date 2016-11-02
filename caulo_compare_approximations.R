library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(scales)

setwd("/home/yannick/simulations/cs_norm")


### read caulobacter dataset and generate with different sampling depths

a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                  skip=0L,nrows=1000000)

#zoom on a portion of the dataset
begin=1000000
end=1500000
begin=2000000
end=2150000
begin=73800287
end=73861120
for (i in c("NcoI","BglIIr1","BglIIr2","BglIIrif")) {
  load(paste0("data/caulo_",i,"_all_csdata_with_data.RData"))
  data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
  cs_data = csnorm:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
  settings=list(circularize=-1,dmin=csd@settings$dmin,dmax=cs_data$biases[,max(pos)-min(pos)]+0.01)
  csd = new("CSdata", info=csd@info, settings=settings,
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  csd@data=data.table()
  save(csd, file=paste0("data/caulo_",i,"_500k_csdata.RData"))
  csd2=csd
}


#normalize with serial sampler
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 1, bf_per_decade = 5, lambdas = 10**seq(from=-1,to=1,length.out=6), ncores = 30, iter = 100000)
cs=run_serial(cs, bf_per_kb = 1, bf_per_decade = 5, init=cs@par, iter = 100000)
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_exact_initgauss.RData")

prefix="tmp/gauss_2enz_500k_lambda"
#prefix="data/caulo_rif_500k_csnorm_optimized_gauss"
#prefix="data/caulo_NcoI_500k_csnorm_optimized_gibbs_simplified_lambda"
#prefix="data/caulo_NcoI_500k_csnorm_optimized_gibbs_gauss_lambda"
registerDoParallel(cores=10)
info=foreach (lambda=10**seq(from=-2,to=2,length.out=10),.combine=rbind, .errorhandling='remove') %dopar% {
  load(paste0(prefix,lambda,".RData"))
  data.table(lambda=lambda,disp.own=cs@par$alpha,val.own=cs@par$value)
}
info[order(val.own)]

#normalize with gibbs sampler
load("data/caulo_NcoI_all_csdata.RData")
csd1=csd
load("data/caulo_BglIIr1_500k_csdata.RData")
csd2=csd
load("data/caulo_BglIIr2_500k_csdata.RData")
csd3=csd
load("data/caulo_BglIIrif_500k_csdata.RData")
csd4=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3,csd4), different.decays=c("enzyme","condition"))
cs=merge_cs_norm_datasets(list(csd1), different.decays="none")
cs = run_gauss_gibbs(cs, bf_per_kb=1, bf_per_decade=20, bins_per_bf=10, init=lambda,
               ngibbs = 1, iter=10000)
cs = run_gauss(cs, bf_per_kb=1, bf_per_decade=10, bins_per_bf=10, ngibbs = 20, iter=100000, init_alpha=1e-7, ncounts = 1000000)
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized.RData")

load("data/caulo_BglII_rif_all_csnorm_optimized.RData")



ggplot(cs@par$biases)+geom_line(aes(pos,log_nu))+
  geom_pointrange(aes(pos,log_nu+z,ymin=log_nu+z+std,ymax=log_nu+z-std,colour=cat), alpha=0.1)+
  facet_wrap(~name)+xlim(1.2e6,1.5e6)+ylim(-2,2)
ggplot(cs@par$biases)+geom_line(aes(pos,log_delta))+
  geom_pointrange(aes(pos,log_delta+z,ymin=log_delta+z+std,ymax=log_delta+z-std,colour=cat), alpha=0.1)+
  facet_wrap(~name)+xlim(1.2e6,1.5e6)

ggplot(cs@par$decay)+geom_line(aes(dist,log_decay),colour="red")+ facet_wrap(~name) +
  geom_pointrange(aes(x=dist,y=log_decay+z,ymin=log_decay+z-std,ymax=log_decay+z+std),alpha=0.1)+
  scale_x_log10()
ggplot(cs@par$decay[dist<2000])+geom_line(aes(dist,log_decay),colour="red")+ facet_wrap(~name) +
  geom_pointrange(aes(x=dist,y=log_decay+z,ymin=log_decay+z-std,ymax=log_decay+z+std),alpha=0.1)
scale_x_log10()




cs=bin_all_datasets(cs, resolution=20000, ncores=30, verbose=T, ice=1)
mat=get_matrices(cs, resolution=20000, group="all")
mat=get_matrices(cs, resolution=20000, group="condition")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(ice.1)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=-log(normalized)))+
  scale_fill_gradient2()+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)
ggsave(filename = "images/caulo_BglII_decay_all_normalized.pdf", width=20, height=20)

a=mat[,.(avg=mean(normalized)),keyby=begin2-begin1]
setkey(a,begin2)
ggplot(a)+geom_point(aes(begin2,avg))



load("data/caulo_NcoI_all_csnorm_optimized.RData")
cs=detect_interactions(cs, resolution=20000, group="all", threshold=0.95, ncores=10)
save(cs, file="data/caulo_NcoI_all_csnorm_optimized.RData")
mat=get_interactions(cs, type="interactions", resolution=20000, group="all", threshold=0.95, ref="expected")
mat[,ncounts:=NULL]
gmat=get_interactions(cs, type="interactions", resolution=20000, group="condition", threshold=0.95, ref="expected")
mat=rbind(mat[name!="rif BglII 1"],gmat[name=="WT"])
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  #geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name, nrow=2, ncol=2)
ggsave(filename = "images/caulo_BglII_decay_all_normalized_interactions.pdf", width=20, height=20)


dmat=get_interactions(cs, type="differences", resolution=20000, group="all",
                      threshold=0.95, ref="WT BglII 2")
dmat[,ncounts:=NULL]
setnames(dmat, "prob.gt.WT BglII 2", "prob.gt.WT")
dgmat=get_interactions(cs, type="differences", resolution=20000, group="condition", threshold=0.95, ref="WT")
dgmat=rbind(dmat[name=="rif BglII 1"],dgmat)
ggplot(dgmat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  geom_point(aes(begin2,begin1,colour=direction),data=dgmat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  facet_grid(~name)+theme_bw()+theme(legend.position = "none")


ggplot(dmat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=-log(signal)), data=dmat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_gradient2()+
  facet_grid(~name)+theme_bw()+theme(legend.position = "none")
ggsave(filename = "images/caulo_BglII_decay_all_differences.pdf", width=22, height=8)





cts=rbind(ctspred[,.(name,id1,id2,pos1,pos2,distance,log_decay,count=contact.close,mean=exp(log_mean_cclose),cat="close")],
             ctspred[,.(name,id1,id2,pos1,pos2,distance,log_decay,count=contact.far,mean=exp(log_mean_cfar),cat="far")],
             ctspred[,.(name,id1,id2,pos1,pos2,distance,log_decay,count=contact.up,mean=exp(log_mean_cup),cat="up")],
             ctspred[,.(name,id1,id2,pos1,pos2,distance,log_decay,count=contact.down,mean=exp(log_mean_cdown),cat="down")])
cts[,bin:=list(cut(distance, begins, ordered_result=T, right=F, include.lowest=T,dig.lab=12))]
ggplot(cts)+geom_jitter(aes(bin,count/mean))

mat[bin1==min(bin1)&bin2==max(bin2)]
cts[bin1==mbin1&bin2==mbin2,.(observed=sum(count),expected=sum(mean))]
cts[bin1==mbin1&bin2==mbin2&mean>1e10]


load("data/caulo_2enz_500k_csnorm_optimized_gauss.RData")
cs = run_serial(cs, init=cs@par, bf_per_kb=1, bf_per_decade=20,  iter=1000, init_alpha=1e-8)
save(cs,file="data/caulo_NcoI_500k_csnorm_optimized_exact_initgauss.RData")
cs = run_gauss_gibbs(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=100, init=cs@par,
               ngibbs = 10, iter=10000)
save(cs,file="data/caulo_2enz_500k_csnorm_optimized_gauss_initexact.RData")




#plots
dsets=c("data/caulo_NcoI_500k_csnorm_optimized_exact_initgauss.RData",
        "data/caulo_NcoI_500k_csnorm_optimized.RData")
names=c("exact",
        "approximation")

#iota and rho
iota = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(name=cs@biases[,name],pos=cs@biases[,pos],iota=exp(cs@par$log_iota),rho=exp(cs@par$log_rho),method=j)
}
ggplot(iota)+geom_line(aes(pos,iota,colour=method))#+facet_wrap(~name)
ggsave(filename = "images/caulo_NcoI_500k_iota_bias.pdf", width=10, height=7)
ggplot(iota)+geom_line(aes(pos,rho,colour=method))#+facet_wrap(~name)
ggsave(filename = "images/caulo_NcoI_500k_rho_bias.pdf", width=10, height=7)
#
ggplot(merge(iota[method=="exact",.(pos,iotaref=iota,rhoref=rho)],iota[method!="exact"],by="pos"))+
  geom_point(aes(iotaref,iota),colour="red")+stat_function(fun=identity)+xlab("exact")+ylab("approximation")
ggsave(filename = "images/caulo_NcoI_500k_iota_bias_correlation.pdf", width=10, height=7)
ggplot(merge(iota[method=="exact",.(pos,iotaref=iota,rhoref=rho)],iota[method!="exact"],by="pos"))+
  geom_point(aes(rhoref,rho),colour="red")+stat_function(fun=identity)+xlab("exact")+ylab("approximation")
ggsave(filename = "images/caulo_NcoI_500k_rho_bias_correlation.pdf", width=10, height=7)
#
cor.test(iota[method=="exact",log(iota)],iota[method=="simplified",log(iota)])
cor.test(iota[method=="exact",log(iota)],iota[method=="approximation",log(iota)])
cor.test(iota[method=="exact",log(rho)],iota[method=="simplified",log(rho)])
cor.test(iota[method=="exact",log(rho)],iota[method=="approximation",log(rho)])

#decay
decay = rbind({load(dsets[1]); cs@par$decay[,.(method="exact",name,dist,decay)]},
              {load(dsets[2]); cs@par$decay[,.(method="approximation",name,dist=distance,decay=exp(kappa-cs@par$eC))]})
ggplot(decay)+geom_line(aes(dist,decay,colour=method))+scale_x_log10()+scale_y_log10()
ggsave(filename = "images/caulo_NcoI_500k_diagonal_decay.pdf", width=10, height=7)
#
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(dist=cs@counts[,distance], log_decay=cs@par$log_decay, method=j)
}
cor.test(decay[method=="exact",log_decay],decay[method=="simplified",log_decay])
cor.test(decay[method=="exact",log_decay],decay[method=="gauss",log_decay])


#parameters
params = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(method=j,name=cs@design[,name],eC=cs@par$eC,alpha=cs@par$alpha,lambda_nu=cs@par$lambda_nu,
             lambda_delta=cs@par$lambda_delta,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
}
setkey(params, name,method)
params

