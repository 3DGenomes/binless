library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


a=examine_dataset("/scratch/ralph/HiC/3_Mapped/EScell_Rbfox1_3.5M.tsv")
csd1=read_and_prepare("/scratch/ralph/HiC/3_Mapped/EScell_Rbfox1_3.5M.tsv",
                      "data/ralph_EScell_Rbfox1_3.5M", "WT", "1", enzyme="MboI",
                      circularize=-1, dangling.L=c(4),
                      dangling.R=c(-1), maxlen=500, save.data=T)


load("data/ralph_EScell_Rbfox1_3.5M_csdata_with_data.RData")
data=csd@data
plot_binned(data, resolution=10000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])
plot_raw(data, b1=1081000, e1=1099000)

plot_binned(data, resolution=1000, b1=1089000, e1=1092000)
plot_binned(data, resolution=100, b1=1089000, e1=1092000)
plot_raw(data, b1=1089000, e1=1094000)

binned = csnorm:::bin_data(csd@data, resolution, b1=csd@biases[,min(pos)], b2=csd@biases[,min(pos)])
ggplot(binned[bin1>300&bin2>300])+geom_raster(aes(bin1,bin2, fill=log(N)))#+
  geom_raster(aes(begin2,begin1, fill=log(N)))+
  xlim(binned[,min(begin1)],binned[,max(begin1)])+
  ylim(binned[,min(begin1)],binned[,max(begin1)])#+
  theme_bw()+
  scale_fill_gradient(low="white", high="black")



### subsample
begin=8500000 #450k
begin=8854480 #100k
end=8954480
load("data/ralph_EScell_Rbfox1_3.5M_csdata_with_data.RData")
data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
cs_data = csnorm:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
csd = new("CSdata", info=csd@info, settings=list(circularize=-1),
          data=data, biases=cs_data$biases, counts=cs_data$counts)
csd@data=data.table()
save(csd, file="data/ralph_EScell_Rbfox1_100k_csdata.RData")


#normalize with serial sampler
load("data/ralph_EScell_Rbfox1_60k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 5, bf_per_decade = 5, lambdas = 10**seq(from=-1,to=1,length.out=6), ncores = 30, iter = 100000)
cs=run_serial(cs, bf_per_kb = 5, bf_per_decade = 5, init=lambda, iter = 1000000)
save(cs, file="data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact.RData")

#normalize with gibbs sampler
load("data/ralph_EScell_Rbfox1_20k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_simplified(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, groups=10, lambdas=10**seq(from=-3,to=0,length.out=8),
                    ngibbs = 20, iter=10000, ncores=30)
cs = run_simplified_gibbs(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, groups=10, init=cs@par,
                          ngibbs = 20, iter=10000)
save(cs, file="data/ralph_EScell_Rbfox1_20k_csnorm_optimized_gibbs_simplified_initexact.RData")

#normalize with gauss sampler
load("data/ralph_EScell_Rbfox1_60k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, lambdas=10**seq(from=-2,to=2,length.out=10),
               ngibbs = 20, iter=10000, ncores=30)
cs = run_gauss_gibbs(cs, bf_per_kb=5, bf_per_decade=5, bins_per_bf=100, init=cs@par, fit.disp=T,
                     ngibbs = 5, iter=10000)
save(cs, file="data/ralph_EScell_Rbfox1_60k_csnorm_optimized_gibbs_gauss_redo.RData")

load("data/rao_HiCall_chrX_60k_csnorm_optimized_gibbs_gauss.RData")
cs=run_serial(cs, bf_per_kb = 5, bf_per_decade = 5, init=cs@par, iter = 10)
save(cs, file="rao_HiCall_chrX_60k_csnorm_optimized_exact_initgauss.RData")


prefix="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_exact_lambda"
#prefix="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_simplified_lambda"
#prefix="data/ralph_EScell_Rbfox1_100k_csnorm_optimized_gauss_lambda"
registerDoParallel(cores=15)
info=foreach (lambda=10**seq(from=-3,to=3,length.out=15),.combine=rbind, .errorhandling='remove') %dopar% {
  load(paste0(prefix,lambda,"_sub10.RData"))
  #out=cs@diagnostics$out
  out=cs@diagnostics$out.init
  #out=cs@diagnostics$out.bias20
  #out=cs@diagnostics$out.decay20
  #out=cs@diagnostics$out.disp20
  #data.table(lambda=lambda,lambda.nu=cs@par$lambda_nu,disp.own=cs@par$alpha,val.own=cs@par$value,out=tail(out,1))
}
info[order(val.own)]

prefix="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_exact_lambda"
#prefix="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_simplified_lambda"
#prefix="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_gauss_lambda"
registerDoParallel(cores=30)
info=foreach (i=Sys.glob(paste0(prefix,"*.RData")),.combine=rbind, .errorhandling='remove') %dopar% {
  load(i)
  out=cs@diagnostics$out
  #out=cs@diagnostics$out.bias20
  #out=cs@diagnostics$out.decay20
  #out=cs@diagnostics$out.disp20
  data.table(fname=i,lambda.nu=cs@par$lambda_nu,disp.own=cs@par$alpha,val.own=cs@par$value,out=tail(out,1))
}
info[order(val.own)]


#plots
dsets=c("data/ralph_EScell_Rbfox1_450k_csnorm_optimized_exact.RData",
        "data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_simplified.RData",
        "data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_gauss.RData")
names=c("exact",
        "simplified",
        "gauss")

dsets=c("data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact.RData",
        "data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact_extension.RData",
        "data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact_initgauss.RData",
        "data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact_initsimplified.RData")
names=c("exact",
        "exacte",
        "exactig",
        "exactis")



#nu and delta
nu = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],nu=exp(cs@par$log_nu),delta=exp(cs@par$log_delta),method=j)
}
ggplot(nu)+geom_line(aes(pos,nu,colour=method))+geom_point(aes(pos,nu),colour="red",data=nu[method=="exact"])
ggplot(nu)+geom_line(aes(pos,delta,colour=method))+geom_point(aes(pos,delta),colour="red",data=nu[method=="exact"])
#
ggplot(merge(nu[method=="exact",.(pos,nuref=nu,deltaref=delta)],nu[method!="exact"],by="pos"))+
  geom_point(aes(nuref,nu,colour=method))+stat_function(fun=identity)
ggplot(merge(nu[method=="exact",.(pos,nuref=nu,deltaref=delta)],nu[method!="exact"],by="pos"))+
  geom_point(aes(deltaref,delta,colour=method))+stat_function(fun=identity)

#nu with points
load(dsets[3])
ggplot(cs@par$biases[id>130])+geom_line(aes(pos,log_nu))+
  geom_pointrange(aes(x=pos,y=log_nu+z,ymin=log_nu+z-std,ymax=log_nu+z+std,colour=cat))

#beta nu
a=data.table(bnu=cs@par$beta_nu[1,])
a[,id:=.I]
a[,pos:=seq(cs@biases[,min(pos)],cs@biases[,max(pos)],length.out=.N)]
ggplot(a)+geom_point(aes(id,bnu))
ggplot(nu)+geom_line(aes(pos,log(nu),colour=method))+geom_line(aes(pos,bnu),data=a)

#logp with gibbs step
load(dsets[2])
logp = foreach(i=0:20, .combine=rbind) %do% {
  if (i==0) {
    data.table(step=i,leg="bias",value=cs@diagnostics$op.init$value)
  } else {
    a=foreach (leg=c("decay","bias","disp"),.combine=rbind) %do% {
      data.table(step=i,leg=leg,value=cs@diagnostics[[paste0("op.",leg,i)]]$value)
    }
    rbind(a,data.table(step=i,leg="lnu",value=log(cs@diagnostics[[paste0("op.",leg,i)]]$par$lambda_nu)))
  }
}
ggplot(logp)+geom_line(aes(step,value))+facet_wrap(~leg, scales="free")

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
  if ("decay" %in% names(cs@par$decay)) {
    cs@par$decay[,.(method=j,dist,decay)]
  } else {
    cs@par$decay[,.(method=j,dist,decay=exp(log_decay))]
  }
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
dset=dsets[2]
ref=dsets[1]
nsteps=20
load(ref)
refparams=data.table(method="ref",step=0:nsteps,leg="ref",eC=cs@par$eC,alpha=cs@par$alpha,
                     lambda_nu=cs@par$lambda_nu,log_nu38=cs@par$log_nu[38],
                     lambda_delta=cs@par$lambda_delta,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
dset="data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact_lambda0.001.RData"
dset="data/ralph_EScell_Rbfox1_100k_csnorm_optimized_gibbs_simplified_lambda100.RData"
dset="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_gauss_lambda100.RData"
load(dset)
params = foreach(i=0:nsteps,.combine=rbind) %do% {
                 if (i==0) {
                   params=data.table(method="10",step=0,leg="init",eC=cs@diagnostics$op.init$par$eC, alpha=cs@diagnostics$op.init$par$alpha,
                                     lambda_nu=cs@diagnostics$op.init$par$lambda_nu, log_nu38=cs@diagnostics$op.init$par$log_nu[38],
                                     beta_diag5=cs@diagnostics$op.init$par$beta_diag[5], lambda_delta=cs@diagnostics$op.init$par$lambda_delta,
                                     lambda_diag=cs@diagnostics$op.init$par$lambda_diag, value=cs@diagnostics$op.init$value, out=tail(cs@diagnostics$out.init,1))
                 } else {
                   foreach (leg=c("decay","bias","disp"), inc=(c(0,0.3,0.6)),.combine=function(...){rbind(...,fill=T)}) %do% {
                     stepname=paste0(leg,i)
                     op=cs@diagnostics[[paste0("op.",stepname)]]
                     as.data.table(c(list(method="10",step=i+inc,leg=leg, log_nu38=op$par$log_nu[38], beta_diag5=op$par$beta_diag[5]),
                                        op$par[c("eC","alpha","lambda_nu","lambda_delta", "lambda_diag")],
                                        list(value=op$value, out=tail(cs@diagnostics[[paste0("out.",stepname)]],1))))
                   }
                 }
}

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
ggplot(params[!(leg%in%c("ref","init"))])+geom_line(aes(step,value))+
  geom_point(aes(step,value,colour=out))+facet_wrap(~leg, scales = "free")
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

for (i in 1:20) {
  cat(i,"decay",tail(cs@diagnostics[[paste0("out.decay",i)]],1),"\n")
  #cat(i,"bias",tail(cs@diagnostics[[paste0("out.bias",i)]],1),"\n")
  #cat(i,"disp",tail(cs@diagnostics[[paste0("out.disp",i)]],1),"\n")
}

#prefix="data/ralph_EScell_Rbfox1_60k_csnorm_optimized_exact_lambda"
#prefix="data/ralph_EScell_Rbfox1_60k_csnorm_optimized_gibbs_simplified_lambda"
prefix="data/ralph_EScell_Rbfox1_450k_csnorm_optimized_gibbs_gauss_lambda"
registerDoParallel(cores=10)
params=foreach (lambda=10**seq(from=-2,to=2,length.out=10),.combine=rbind, .errorhandling='remove') %dopar% {
  load(paste0(prefix,lambda,".RData"))
  params = foreach(i=0:nsteps,.combine=rbind) %do% {
    if (i==0) {
      params=data.table(method="10",step=0,leg="init",runtime=cs@diagnostics$runtime.init,eC=cs@diagnostics$op.init$par$eC, alpha=cs@diagnostics$op.init$par$alpha,
                        lambda_nu=cs@diagnostics$op.init$par$lambda_nu, log_nu38=cs@diagnostics$op.init$par$log_nu[38],
                        beta_diag5=cs@diagnostics$op.init$par$beta_diag[5], lambda_delta=cs@diagnostics$op.init$par$lambda_delta,
                        lambda_diag=cs@diagnostics$op.init$par$lambda_diag, value=cs@diagnostics$op.init$value, out=tail(cs@diagnostics$out.init,1))
    } else {
      foreach (leg=c("decay","bias","disp"), inc=(c(0,0.3,0.6)),.combine=function(...){rbind(...,fill=T)}) %do% {
        stepname=paste0(leg,i)
        op=cs@diagnostics[[paste0("op.",stepname)]]
        as.data.table(c(list(method="10",step=i+inc,leg=leg,log_nu38=op$par$log_nu[38], beta_diag5=op$par$beta_diag[5]),
                        op$par[c("eC","alpha","lambda_nu","lambda_delta", "lambda_diag")],
                        list(value=op$value, out=tail(cs@diagnostics[[paste0("out.",stepname)]],1),
                             runtime=cs@diagnostics[[paste0("runtime.",stepname)]])))
      }
    }
  }
  params[,lambda:=lambda]
  params
}
#value
ggplot(params[!(leg%in%c("ref","init"))])+geom_line(aes(step,value))+
  geom_point(aes(step,value,colour=out))+facet_grid(leg~lambda, scales = "free")
#dispersion
ggplot(params[!(leg%in%c("ref","init"))])+geom_line(aes(step,value))+
  geom_point(aes(step,alpha,colour=out))+facet_grid(leg~lambda, scales = "free")
#lambda_nu
ggplot(params[!(leg%in%c("ref","init"))])+geom_line(aes(step,value))+
  geom_point(aes(step,alpha,colour=out))+facet_grid(leg~lambda, scales = "free")
#runtime
ggplot(params[!(leg%in%c("ref","init"))])+geom_line(aes(step,runtime))+
  geom_point(aes(step,runtime,colour=out))+facet_grid(leg~lambda, scales = "free")


