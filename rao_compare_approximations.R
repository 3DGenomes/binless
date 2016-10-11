library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


### read caulobacter dataset and generate with different sampling depths


#normalize with serial sampler
load("data/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 1, bf_per_decade = 5, lambdas = 10**seq(from=-1,to=1,length.out=6), ncores = 30, iter = 100000)
cs=run_serial(cs, bf_per_kb = 0.25, bf_per_decade = 5, init=cs@par, iter = 100000)
save(cs, file="data/rao_HiCall_chrX_450k_csnorm_optimized_exact.RData")

for (lambda in c(round(10**seq(from=-3,to=0,length.out=8),digits=3),3.000,7.000)) {
  try(load(paste0("data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_gauss_",lambda,".RData")))
  cat(lambda," ",cs@par$value," ",cs@par$alpha,"\n")
}

#normalize with gibbs sampler
load("data/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, groups=10, lambdas=10**seq(from=-3,to=0,length.out=8),
                    ngibbs = 20, iter=10000, ncores=30)
cs = run_simplified_gibbs(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, groups=10, init=cs@par,
                          ngibbs = 20, iter=10000)
save(cs, file="data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group10_initexact.RData")

#normalize with gauss sampler
load("data/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, lambdas=10**seq(from=-3,to=0,length.out=8),
               ngibbs = 3, iter=10000, ncores=30)
cs = run_gauss_gibbs(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, init=1,
                     ngibbs = 5, iter=10000)
save(cs, file="data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_gauss_initexact.RData")

nu=rbind(data.table(pos=cs@biases[,pos],method="1",nu=exp(op.1$par$log_nu),delta=exp(op.1$par$log_delta)),
         data.table(pos=cs@biases[,pos],method="10",nu=exp(op.10$par$log_nu),delta=exp(op.10$par$log_delta)),
         data.table(pos=cs@biases[,pos],method="gauss",nu=exp(op.gauss$par$log_nu),delta=exp(op.gauss$par$log_delta)))
ggplot(nu)+geom_line(aes(pos,nu,colour=method))
ggplot(nu)+geom_line(aes(pos,delta,colour=method))

#plots
dsets=c("data/rao_HiCall_chrX_450k_csnorm_optimized_exact.RData",
        #"data/rao_HiCall_chrX_450k_csnorm_optimized_exact_initgauss.RData",
        "data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group10.RData",
        #"data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group10_initexact.RData",
        #"data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group1.RData",
        #"data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group1_initexact.RData",
        "data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_gauss.RData",
        "data/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_gauss_initexact.RData")

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


load(dsets[3])
ggplot(cs@par$biases[id<100])+geom_line(aes(pos,log_nu),colour="red")+
  geom_pointrange(aes(x=pos,y=log_nu+z,ymin=log_nu+z-std,ymax=log_nu+z+std,colour=cat))
ggplot(cs@par$biases)+geom_line(aes(pos,log_delta),colour="red")+
  geom_pointrange(aes(x=pos,y=log_delta+z,ymin=log_delta+z-std,ymax=log_delta+z+std,colour=cat),alpha=0.1)+
  scale_x_log10()


#nu and delta with gibbs step
load(dsets[3])
nu = foreach(i=0:20, .combine=rbind) %do% {
  if (i==0) op=cs@diagnostics$op.init else op=cs@diagnostics[[paste0("op.bias",i)]]
  data.table(pos=cs@biases[,pos],round=i,nu=exp(op$par$log_nu),delta=exp(op$par$log_delta))
}
ggplot(nu)+geom_line(aes(pos,nu))+facet_grid(round~.)
ggplot(nu)+geom_line(aes(pos,delta))+facet_grid(round~.)

#decay
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  cs@par$decay[,.(method=j,dist,decay)]
}
ggplot(decay)+geom_line(aes(dist,decay,colour=method))+scale_x_log10()+scale_y_log10()
ggplot(decay)+geom_line(aes(dist,decay*dist^0.5,colour=method))+scale_x_log10()+scale_y_log10()

load(dsets[3])
ggplot(cs@par$decay)+geom_line(aes(dist,log_decay),colour="red")+
  geom_pointrange(aes(x=dist,y=log_decay+z,ymin=log_decay+z-std,ymax=log_decay+z+std),alpha=0.1)+
  scale_x_log10()


#parameters
params = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(method=j,eC=cs@par$eC,alpha=cs@par$alpha,lambda_nu=cs@par$lambda_nu,
             lambda_delta=cs@par$lambda_delta,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
}
params

#parameters with gibbs step
dset=dsets[3]
ref=dsets[3]
nsteps=20
load(ref)
refparams=data.table(method="ref",step=0:nsteps,leg="ref",eC=cs@par$eC,alpha=cs@par$alpha,
                     lambda_nu=cs@par$lambda_nu,log_nu38=cs@par$log_nu[38],
                     lambda_delta=cs@par$lambda_delta,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
load(dset)
params = rbind(params,
               foreach(i=0:nsteps,.combine=rbind) %do% {
                 if (i==0) {
                   params=data.table(method="10",step=0,leg="init",eC=cs@diagnostics$op.init$par$eC, alpha=cs@diagnostics$op.init$par$alpha,
                                     lambda_nu=cs@diagnostics$op.init$par$lambda_nu, log_nu38=cs@diagnostics$op.init$par$log_nu[38],
                                     beta_diag5=cs@diagnostics$op.init$par$beta_diag[5], lambda_delta=cs@diagnostics$op.init$par$lambda_delta,
                                     lambda_diag=cs@diagnostics$op.init$par$lambda_diag, value=cs@diagnostics$op.init$value)
                 } else {
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
                 }
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
ggplot(params)+geom_line(aes(step,value))+geom_point(aes(step,value,colour=leg))+
  geom_hline(aes(yintercept=value),data=refparams,colour="red")
#alpha
ggplot(params)+geom_line(aes(step,alpha))+geom_point(aes(step,alpha,colour=leg))+
  geom_hline(aes(yintercept=alpha),data=refparams,colour="red")

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
load("data/old/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 1, bf_per_decade = 5, lambdas = 10**seq(from=-2,to=0,length.out=6), ncores = 30, iter = 100000)
save(cs, file="data/old/rao_HiCall_chrX_450k_csnorm_optimized_exact.RData")

#normalize with gibbs sampler
load("data/old/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=10, lambdas=10**seq(from=-2,to=0,length.out=6),
                    ngibbs = 5, iter=10000, ncores=30)
save(cs, file="data/old/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group10.RData")

#normalize with gibbs sampler
load("data/old/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=1, bf_per_decade=5, bins_per_bf=10, groups=1, lambdas=10**seq(from=-2,to=0,length.out=6),
                    ngibbs = 5, iter=10000, ncores=30)
save(cs, file="data/old/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group1.RData")




#plots
dsets=c("data/old/rao_HiCall_chrX_450k_csnorm_optimized_exact.RData",
        "data/old/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group10.RData",
        "data/old/rao_HiCall_chrX_450k_csnorm_optimized_gibbs_quantile_group1.RData")
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


