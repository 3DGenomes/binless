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
begin=500000
end=1000000
begin=2000000
end=2150000
begin=73800287
end=73861120
for (i in c("BglIIr1","BglIIr2","BglII_rifampicin")) {
  load(paste0("data/caulo_",i,"_all_csdata_with_data.RData"))
  data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
  cs_data = csnorm:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
  csd = new("CSdata", info=csd@info, settings=list(circularize=-1),
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  save(csd, file=paste0("data/caulo_",i,"_500k_csdata_with_data.RData"))
  csd@data=data.table()
  save(csd, file=paste0("data/caulo_",i,"_500k_csdata.RData"))
  csd2=csd
}


#normalize with serial sampler
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 1, bf_per_decade = 5, lambdas = 10**seq(from=-1,to=1,length.out=6), ncores = 30, iter = 100000)
cs=run_serial(cs, bf_per_kb = 0.25, bf_per_decade = 5, init=cs@par, iter = 100000)
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_exact.RData")

for (lambda in c(0.01,0.05,0.1,0.5,1,5,10)) {
  try(load(paste0("data/caulo_NcoI_500k_csnorm_optimized_exact_lambda",lambda,".RData")))
  cat(lambda," ",cs@par$value,"\n")
}


#normalize with gibbs sampler
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=-1,to=1,length.out=6),
                    ngibbs = 20, iter=10000, ncores=30)
cs = run_simplified_gibbs(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, init=refpar,
                    ngibbs = 20, iter=10000)
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_gibbs_quantile_group10.RData")

#normalize with gibbs sampler
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=1, lambdas=10**seq(from=-1,to=1,length.out=6),
                    ngibbs = 20, iter=10000, ncores=30)
cs = run_simplified_gibbs(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=1, init=refpar,
                          ngibbs = 20, iter=10000)
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_gibbs_quantile_group1.RData")

#normalize with gauss sampler
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=100, lambdas=10**seq(from=-1,to=1,length.out=6),
                    ngibbs = 20, iter=10000, ncores=30)
cs = run_gauss_gibbs(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=100, init=cs@par,
               ngibbs = 20, iter=10000)
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_gibbs_gauss_initexact.RData")


nu=rbind(data.table(pos=cs@biases[,pos],method="1",nu=exp(op.1$par$log_nu),delta=exp(op.1$par$log_delta)),
         data.table(pos=cs@biases[,pos],method="10",nu=exp(op.10$par$log_nu),delta=exp(op.10$par$log_delta)),
         data.table(pos=cs@biases[,pos],method="gauss",nu=exp(op.gauss$par$log_nu),delta=exp(op.gauss$par$log_delta)))
ggplot(nu)+geom_line(aes(pos,nu,colour=method))
ggplot(nu)+geom_line(aes(pos,delta,colour=method))

#plots
dsets=c("data/caulo_NcoI_500k_csnorm_optimized_exact.RData",
        #"data/caulo_NcoI_500k_csnorm_optimized_exact_initgauss.RData",
        "data/caulo_NcoI_500k_csnorm_optimized_gibbs_quantile_group10.RData",
        #"data/caulo_NcoI_500k_csnorm_optimized_gibbs_quantile_group10_initexact.RData",
        #"data/caulo_NcoI_500k_csnorm_optimized_gibbs_quantile_group1.RData",
        #"data/caulo_NcoI_500k_csnorm_optimized_gibbs_quantile_group1_initexact.RData",
        "data/caulo_NcoI_500k_csnorm_optimized_gibbs_gauss.RData",
        "data/caulo_NcoI_500k_csnorm_optimized_gibbs_gauss_initexact.RData")
        
names=c("exact",
        #"exactig",
        "group10",
        #"group10i",
        #"group1",
        #"group1i",
        "gauss",
        "gaussi")
        


#nu and delta
nu = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],nu=exp(cs@par$log_nu),delta=exp(cs@par$log_delta),method=j)
}
ggplot(nu)+geom_line(aes(pos,nu,colour=method))+geom_point(aes(pos,nu),colour="red",data=nu[method=="exact"])
ggplot(nu)+geom_line(aes(pos,delta,colour=method))+geom_point(aes(pos,delta),colour="red",data=nu[method=="exact"])

#decay
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  cs@par$decay[,.(method=j,dist,decay)]
}
ggplot(decay)+geom_line(aes(dist,decay,colour=method))+scale_x_log10()+scale_y_log10()
ggplot(decay)+geom_line(aes(dist,decay*dist^0.5,colour=method))+scale_x_log10()+scale_y_log10()

load(dsets[3])
ggplot(cs@par$decay)+geom_line(aes(dist,decay),colour="red")+
  geom_errorbar(aes(dist,ymax=exp(kappa_hat-cs@par$eC)+sdl,ymin=exp(kappa_hat-cs@par$eC)-sdl),alpha=0.1)+
  geom_point(aes(dist,exp(kappa_hat-cs@par$eC)),alpha=0.1)+
  scale_x_log10()+scale_y_log10()


#parameters
params = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(method=j,eC=cs@par$eC,alpha=cs@par$alpha,lambda_nu=cs@par$lambda_nu,
             lambda_delta=cs@par$lambda_delta,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
}
params

#parameters with gibbs step
dset=dsets[2]
ref=dsets[1]
nsteps=20
load(ref)
refparams=data.table(method="ref",step=0:nsteps,leg="ref",eC=cs@par$eC,alpha=cs@par$alpha,
                     lambda_nu=cs@par$lambda_nu,log_nu38=cs@par$log_nu[38],
                     lambda_delta=cs@par$lambda_delta,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
load(dset)
cs@diagnostics$op.init$value=refparams$value
params=data.table(method="10",step=0,leg="bias",eC=cs@diagnostics$op.init$par$eC, alpha=cs@diagnostics$op.init$par$alpha,
                 lambda_nu=cs@diagnostics$op.init$par$lambda_nu, log_nu38=cs@diagnostics$op.init$par$log_nu[38],
                 beta_diag5=cs@diagnostics$op.init$par$beta_diag[5], lambda_delta=cs@diagnostics$op.init$par$lambda_delta,
                 lambda_diag=cs@diagnostics$op.init$par$lambda_diag, value=cs@diagnostics$op.init$value)
params = rbind(params,
               foreach(i=1:nsteps,.combine=rbind) %do% {
  stepname=paste0("decay",i)
  op=cs@diagnostics[[paste0("op.",stepname)]]
  a1=as.data.table(c(list(method="10",step=i,leg="decay", log_nu38=op$par$log_nu[38], beta_diag5=op$par$beta_diag[5]),
                     op$par[c("eC","alpha","lambda_nu","lambda_delta", "lambda_diag")],
               list(value=op$value)))
  #
  stepname=paste0("bias",i)
  op=cs@diagnostics[[paste0("op.",stepname)]]
  a2=as.data.table(c(list(method="10",step=i+0.3,leg="bias", log_nu38=op$par$log_nu[38], beta_diag5=op$par$beta_diag[5]),
                     op$par[c("eC","alpha","lambda_nu","lambda_delta", "lambda_diag")],
                  list(value=op$value)))
  #
  stepname=paste0("disp",i)
  op=cs@diagnostics[[paste0("op.",stepname)]]
  a3=as.data.table(c(list(method="10",step=i+0.6,leg="disp", log_nu38=op$par$log_nu[38], beta_diag5=op$par$beta_diag[5]),
                     op$par[c("eC","alpha","lambda_nu","lambda_delta", "lambda_diag")],
                  list(value=op$value)))
  rbind(a1,a2,a3,fill=T)
})

#eC
ggplot(params)+geom_line(aes(step,eC))+geom_point(aes(step,eC,colour=leg))+
  geom_hline(aes(yintercept=eC),data=refparams,colour="red")
#lambda_nu
ggplot(params)+geom_line(aes(step,lambda_nu))+geom_point(aes(step,lambda_nu,colour=leg))+
  geom_hline(aes(yintercept=lambda_nu),data=refparams,colour="red")
#log_nu38
ggplot(params)+geom_line(aes(step,log_nu38))+geom_point(aes(step,log_nu38,colour=leg))+
  geom_hline(aes(yintercept=log_nu38),data=refparams,colour="red")
#value
ggplot(params)+geom_line(aes(step,value))+geom_point(aes(step,log_nu38,colour=leg))+
  geom_hline(aes(yintercept=log_nu38),data=refparams,colour="red")

mparams=melt(params,id.vars=c("method","step","leg"))
ggplot(mparams)+geom_line(aes(step,value))+geom_point(aes(step,value,colour=leg))+
  geom_hline(aes(yintercept=value),data=melt(refparams,id.vars=c("method","step","leg")),colour="red")+
  facet_wrap(~variable,scales = "free")

ggplot(mparams[variable%in%c("alpha","log_nu38","lambda_nu")])+ylim(0,1)+
  geom_line(aes(step,value,colour=variable))+geom_point(aes(step,value,colour=variable,shape=leg))


#spline parameters
beta = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(id=1:length(cs@par$beta_nu),beta_nu=cs@par$beta_nu[1,],beta_delta=cs@par$beta_delta[1,],method=j)
}
ggplot(beta)+geom_line(aes(id,beta_nu,colour=method))
ggplot(beta)+geom_line(aes(id,beta_delta,colour=method))
#
beta = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(id=1:length(cs@par$beta_diag),beta_diag=cs@par$beta_diag[1,],method=j)
}
ggplot(beta)+geom_line(aes(id,beta_diag,colour=method))
#




#normalize with serial sampler
load("data/old/rao_HiCall_chrX_60k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 1, bf_per_decade = 5, lambdas = 10**seq(from=-2,to=0,length.out=6), ncores = 30, iter = 100000)
save(cs, file="data/old/rao_HiCall_chrX_60k_csnorm_optimized_exact.RData")

#normalize with gibbs sampler
load("data/old/rao_HiCall_chrX_60k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=-2,to=0,length.out=6),
                    ngibbs = 5, iter=10000, ncores=30)
save(cs, file="data/old/rao_HiCall_chrX_60k_csnorm_optimized_gibbs_quantile_group10.RData")

#normalize with gibbs sampler
load("data/old/rao_HiCall_chrX_60k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=1, lambdas=10**seq(from=-2,to=0,length.out=6),
                    ngibbs = 5, iter=10000, ncores=30)
save(cs, file="data/old/rao_HiCall_chrX_60k_csnorm_optimized_gibbs_quantile_group1.RData")




#plots
dsets=c("data/old/rao_HiCall_chrX_60k_csnorm_optimized_exact.RData",
        "data/old/rao_HiCall_chrX_60k_csnorm_optimized_gibbs_quantile_group10.RData",
        "data/old/rao_HiCall_chrX_60k_csnorm_optimized_gibbs_quantile_group1.RData")
names=c("exact","group10","group1")
nu = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],nu=exp(cs@par$log_nu),delta=exp(cs@par$log_delta),method=j)
}
ggplot(nu)+geom_line(aes(pos,nu,colour=method))
ggplot(nu)+geom_line(aes(pos,delta,colour=method))
#
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  cs@par$decay[,.(method=j,dist,decay)]
}
ggplot(decay)+geom_line(aes(dist,decay,colour=method))+scale_x_log10()+scale_y_log10()


