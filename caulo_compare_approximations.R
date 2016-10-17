library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


### read caulobacter dataset and generate with different sampling depths

a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                  skip=0L,nrows=1000000)

#zoom on a portion of the dataset
begin=1000000
end=1150000
begin=2000000
end=2150000
begin=73800287
end=73861120
for (i in c("NcoI","BglIIr1","BglIIr2","BglIIrif")) {
  load(paste0("data/caulo_",i,"_all_csdata_with_data.RData"))
  data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
  cs_data = csnorm:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
  csd = new("CSdata", info=csd@info, settings=list(circularize=-1),
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
load("data/caulo_NcoI_500k_csdata.RData")
csd1=csd
load("data/caulo_BglIIr1_500k_csdata.RData")
csd2=csd
load("data/caulo_BglIIr2_500k_csdata.RData")
csd3=csd
load("data/caulo_BglIIrif_500k_csdata.RData")
csd4=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3,csd4), different.decays=c("enzyme","condition"))
cs = run_gauss(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=100, lambdas=10**seq(from=-2,to=2,length.out=10),
                    ngibbs = 10, iter=10000, ncores=30, prefix="tmp/gauss_2enz_500k")
save(cs, file="data/caulo_2enz_500k_csnorm_optimized_gauss.RData")


load("data/caulo_2enz_500k_csnorm_optimized_gauss.RData")
cs = run_serial(cs, init=cs@par, bf_per_kb=0.25, bf_per_decade=5,  iter=10000, init_alpha=1e-8)
save(cs,file="data/caulo_2enz_500k_csnorm_optimized_exact_initgauss.RData")
cs = run_gauss_gibbs(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=100, init=cs@par,
               ngibbs = 10, iter=10000)
save(cs,file="data/caulo_2enz_500k_csnorm_optimized_gauss_initexact.RData")




#plots
dsets=c(#"data/caulo_2enz_500k_csnorm_optimized_exact.RData",
        "data/caulo_2enz_500k_csnorm_optimized_exact_initgauss.RData",
        "data/caulo_2enz_500k_csnorm_optimized_gauss_initexact.RData",
        "data/caulo_2enz_500k_csnorm_optimized_gauss.RData")
names=c(#"exact",
        "exactig",
        "gaussie",
        "gauss")

#nu and delta
nu = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(name=cs@biases[,name],pos=cs@biases[,pos],nu=exp(cs@par$log_nu),delta=exp(cs@par$log_delta),method=j)
}
ggplot(nu)+geom_line(aes(pos,nu,colour=method))+facet_wrap(~name)
ggsave(filename = "images/caulo_NcoI_500k_nu_bias.pdf", width=10, height=7)
ggplot(nu)+geom_line(aes(pos,delta,colour=method))+facet_wrap(~name)
ggsave(filename = "images/caulo_NcoI_500k_delta_bias.pdf", width=10, height=7)
#
ggplot(merge(nu[method=="exact",.(pos,nuref=nu,deltaref=delta)],nu[method!="exact"],by="pos"))+
  geom_point(aes(nuref,nu,colour=method))+stat_function(fun=identity)
ggsave(filename = "images/caulo_NcoI_500k_nu_bias_correlation.pdf", width=10, height=7)
ggplot(merge(nu[method=="exact",.(pos,nuref=nu,deltaref=delta)],nu[method!="exact"],by="pos"))+
  geom_point(aes(deltaref,delta,colour=method))+stat_function(fun=identity)
ggsave(filename = "images/caulo_NcoI_500k_delta_bias_correlation.pdf", width=10, height=7)
#
cor.test(nu[method=="exact",log(nu)],nu[method=="simplified",log(nu)])
cor.test(nu[method=="exact",log(nu)],nu[method=="approximation",log(nu)])
cor.test(nu[method=="exact",log(delta)],nu[method=="simplified",log(delta)])
cor.test(nu[method=="exact",log(delta)],nu[method=="approximation",log(delta)])

#decay
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  if ("decay" %in% names(cs@par$decay)) {
    cs@par$decay[,.(method=j,name,dist,decay)]
  } else {
    cs@par$decay[,.(method=j,name,dist,decay=exp(log_decay))]
  }
}
ggplot(decay)+geom_line(aes(dist,decay,colour=method))+facet_grid(method~name)+scale_x_log10()+scale_y_log10()
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

