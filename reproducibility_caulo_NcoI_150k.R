library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(Hmisc)

setwd("/home/yannick/simulations/cs_norm")

### read caulobacter dataset and normalize a 150k square
### - with different bf_per_kb
### - with different initial conditions

#zoom on a portion of the dataset
load("data/caulo_NcoI_all_csdata_with_data.RData")
data=cs@data[re.closest1>=2000000&re.closest1<=2150000&re.closest2>=2000000&re.closest2<=2150000]
cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=4042929)
dset_statistics(cs_data$biases,cs_data$counts)
message("*** WRITE")
csd = new("CSdata", info=cs@info, settings=cs@settings,
         data=data, biases=cs_data$biases, counts=cs_data$counts)
if (save.data==T) save(cs, file="data/caulo_NcoI_150k_csdata_with_data.RData")
cs@data=data.table()
save(cs, file="data/caulo_NcoI_150k_csdata.RData")

bf_per_kb=signif(10^(seq(-1,1,length.out = 9))*1/4,digits=2)
lambdas=c(0.01,0.1,1,10,100)



### reproducibility: influence of bf_per_kb (averaging over lambdas)
registerDoParallel(cores=30)
foreach (bpk=bf_per_kb) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_NcoI_150k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  dmin=1-0.01
  dmax=150000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
                                         counts=cs@counts, biases=cs@biases, lambda=lambda,
                                         bf_per_kb=bpk, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    model=csnorm:::stanmodels$fit, biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
    bf_per_kb=bpk, bf_per_decade=bf_per_decade, iter=100000, verbose = F, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  op$par$value=op$value
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bpk, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_150k_bfpkb",bpk,"_lambda",lambda,"_csnorm_optimized.RData"))
}


#generate plots
prefix="NcoI_150k"

#outputs and runtime
outputs = foreach (i=bf_per_kb,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
  data.table(dset=paste("bf",i,"lam",j),bfpkb=i,lambda=j,
             out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime, iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$value)
}
outputs
summary(lm(data=outputs, log(runtime)~log(bfpkb)+log(lambda)))
ggplot(outputs)+geom_point(aes(bfpkb,runtime,colour=log(lambda)))+scale_x_log10()+scale_y_log10()+
  ylab("run time (s)")+xlab("basis functions per kb")+ggtitle("t ~ bf^0.56")
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_runtime.png"), width=10, height=7.5)




#nu and delta: init
registerDoParallel(cores=30)
nu.init = foreach (i=bf_per_kb,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("bf",i,"lam",j),bfpkb=i,lambda=j,logp=cs@par$value)]
}
nu.init[,rank:=frank(-logp,ties.method="dense"),by=bfpkb]
nu.init=nu.init[rank==1]
#nu.init[,c("bfpkb","lambda"):=list(factor(bfpkb),factor(lambda))]
ggplot(nu.init)+geom_line(aes(pos,exp(log_nu)))+facet_wrap(~bfpkb)+ylim(0,2)

#nu and delta: plots and correlation
registerDoParallel(cores=30)
nu = foreach (i=bf_per_kb,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("bf",i,"lam",j),bfpkb=i,lambda=j,logp=cs@par$value)]
}
nu[,rank:=frank(-logp,ties.method="dense"),by=bfpkb]
nu=nu[rank==1]
#nu[,c("bfpkb","lambda"):=list(factor(bfpkb),factor(lambda))]
#
ggplot(nu)+geom_line(aes(pos,exp(log_nu)))+facet_wrap(~bfpkb)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_nu.png"), width=10, height=7.5)
#nu*delta
ggplot(nu[,.(dset,pos,delta=exp(log_nu+log_delta))])+
  geom_line(aes(pos,delta),colour="black")+facet_wrap(~dset)+ylim(0,2)
#delta
ggplot(nu[,.(dset,pos,delta=exp(log_delta))])+
  geom_line(aes(pos,delta),colour="black")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_delta.png"), width=10, height=7.5)
#fpt
get_fpt=function(x,y){a=acf(y,plot=F,lag.max=10000); l=a$lag[a$acf<=0][1]; x[l+1]-x[1]}
fpts=nu[,.(fpt=get_fpt(pos,exp(log_nu+log_delta))),by=dset]
ggplot(fpts)+geom_point(aes(4/dset,fpt/1000))+
  ylab("first passage time of nu x delta (in kb)")+xlab("span of one basis function (in kb)")+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_nu_delta_fpt.png"), width=10, height=7.5)

#lFC
lFC = foreach (i=bf_per_kb,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=i,lFC)]
}
lFC[,dset:=ordered(dset,levels=bf_per_kb)]
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=bf_per_kb,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=i,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
mat[,dset:=ordered(dset,levels=bf_per_kb)]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_normalized.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=bf_per_kb,.combine=rbind,.export="data.table") %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_csnorm_optimized.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay), dset=i,key="dist")
}
decay[,dset:=ordered(dset,levels=bf_per_kb)]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_decay.png"), width=10, height=7.5)




### reproducibility: influence of lambda (for different sampling depths)

lambdas=signif(10^(seq(-2,2,length.out = 9)),digits=2)
samplings=c(5,40,75,400)

registerDoParallel(cores=30)
foreach (lambda=lambdas) %:% foreach (samp=samplings) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",samp,"k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  bf_per_kb=0.25
  dmin=1-0.01
  dmax=150000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, lambda=lambda, verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    model=csnorm:::stanmodels$fit, biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = F, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  op$par$logp=op$value
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_150k_sub",samp,"k_lambda",lambda,"_csnorm_optimized.RData"))
}


#ggplot(data.table(pos=1:36,beta.init=init.op$beta_nu,beta=op$par$beta_nu))+geom_line(aes(pos,beta.init),colour="red")+geom_line(aes(pos,beta),colour="black")

prefix="NcoI_150k"

#outputs and runtime
outputs = foreach (i=lambdas,.combine=rbind) %:% foreach (j=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_sub",j,"k_lambda",i,"_csnorm_optimized.RData"))
  data.table(dset=paste("samp",j,"per",i),samp=j,lambda=i,
             out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime, iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$logp)
}
setkey(outputs,dset)
outputs
ggplot(outputs)+geom_point(aes(dset,logp))
ggplot(outputs)+geom_point(aes(lambda_nu,logp))+scale_x_log10()


#nu and delta: init
registerDoParallel(cores=30)
nu.init = foreach (i=lambdas,.combine=rbind) %:% foreach (j=samplings,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",j,"k_lambda",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("samp",j,"per",i),samp=j,lambda=i,logp=cs@par$logp)]
}
setkey(nu.init,dset)
nu.init[,rank:=frank(-logp,ties.method="dense"),by=samp]
nu.init[,rank:=ordered(rank)]
#levels(nu.init$rank)<-c("1","2","3",rep(NA,nu.init[,nlevels(rank)-3]))
nu.init[,c("samp","lambda"):=list(factor(samp),factor(lambda))]
ggplot(nu.init)+geom_line(aes(pos,exp(log_nu)))+
  facet_grid(lambda~samp, labeller=label_both)+ylim(0,2)
ggsave(filename=paste0("images/caulo_NcoI_150k_initcond_reproducibility_nu_init.png"), width=10, height=7.5)

#nu and delta: plots and correlation
registerDoParallel(cores=30)
nu = foreach (i=lambdas,.combine=rbind) %:% foreach (j=samplings,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",j,"k_lambda",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("samp",j,"per",i),samp=j,lambda=i,
                                        logp=cs@par$logp,lambda_nu=cs@par$lambda_nu)]
}
setkey(nu,dset)
nu[,rank:=frank(-logp,ties.method="dense"),by=samp]
nu[,rank:=ordered(rank)]
#levels(nu$rank)<-c("1","2","3",rep(NA,nu[,nlevels(rank)-3]))
nu[,c("samp","lambda"):=list(factor(samp),factor(lambda))]
#
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=rank))+facet_grid(lambda~samp, labeller=label_both)+ylim(0,2)
ggsave(filename=paste0("images/caulo_NcoI_150k_initcond_reproducibility_nu.png"), width=10, height=7.5)

ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=log(lambda_nu)))+facet_grid(rank~samp, labeller=label_both)+ylim(0,2)
ggsave(filename=paste0("images/caulo_NcoI_150k_initcond_reproducibility_nu2.png"), width=10, height=7.5)






### reproducibility: fitting with a subset of all points

lambdas=signif(10^(seq(-2,2,length.out = 5)),digits=2)
subs=c(1,5,10,25,50,75,90,100)
#subs=c(1,5,50,75,90,100)

registerDoParallel(cores=30)
foreach (sub=subs) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_NcoI_150k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  bf_per_kb=0.25
  dmin=1-0.01
  dmax=150000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, lambda=lambda, verbose=F, iter=10000)))
  counts.sub=cs@counts[sample(.N,round(sub/100*.N))]
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    model=csnorm:::stanmodels$fit, biases=cs@biases, counts = counts.sub, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = F, init=init.op, weight=1)))#weight=cs@counts[,.N]/counts.sub[,.N])))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  op$par$logp=op$value
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  op$par$counts.sub=counts.sub
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_150k_",sub,"pc_lambda",lambda,"_csnorm_optimized.RData"))
}

prefix="NcoI_150k"

#outputs and runtime
outputs = foreach (i=lambdas,.combine=rbind) %:% foreach (j=subs,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_",j,"pc_lambda",i,"_csnorm_optimized.RData"))
  data.table(dset=paste("sub",j,"per",i),sub=j,lambda=i,
             out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime, iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$logp)
}
setkey(outputs,dset)
outputs
summary(lm(log(runtime)~log(sub),data=outputs))
ggplot(outputs)+geom_point(aes(sub,runtime,colour=logp))


#nu and delta: init
registerDoParallel(cores=30)
nu.init = foreach (i=lambdas,.combine=rbind) %:% foreach (j=subs,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_",j,"pc_lambda",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("samp",j,"per",i),
                                        dset=paste("sub",j,"per",i),sub=j,lambda=i,logp=cs@par$logp)]
}
setkey(nu.init,dset)
nu.init[,rank:=frank(-logp,ties.method="dense"),by=sub]
nu.init[,rank:=ordered(rank)]
ggplot(nu.init)+geom_line(aes(pos,exp(log_nu)))+facet_grid(lambda~sub, labeller=label_both)+ylim(0,2)
ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_init.png"), width=10, height=7.5)

#nu: plot unweighted
registerDoParallel(cores=30)
nunowt = foreach (i=lambdas,.combine=rbind) %:% foreach (j=subs,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_",j,"pcnowt_lambda",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("samp",j,"per",i),lambda_nu=cs@par$lambda_nu,
                                        dset=paste("sub",j,"per",i),sub=j,lambda=i,logp=cs@par$logp)]
}
setkey(nunowt,dset)
nunowt[,rank:=frank(-logp,ties.method="dense"),by=sub]
#nu[,rank:=ordered(rank)]
#
ggplot(nunowt)+geom_line(aes(pos,exp(log_nu)))+facet_grid(rank~sub, labeller=label_both)+ylim(0,2)+ylab("nu")
ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_noweight.png"), width=10, height=7.5)

#nu: plot weighted
registerDoParallel(cores=30)
nu = foreach (i=lambdas,.combine=rbind) %:% foreach (j=subs,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_",j,"pc_lambda",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("samp",j,"per",i),lambda_nu=cs@par$lambda_nu,
                                        dset=paste("sub",j,"per",i),sub=j,lambda=i,logp=cs@par$logp)]
}
setkey(nu,dset)
nu[,rank:=frank(-logp,ties.method="dense"),by=sub]
#nu[,rank:=ordered(rank)]
#
ggplot(nu)+geom_line(aes(pos,exp(log_nu)))+facet_grid(rank~sub, labeller=label_both)+ylim(0,2)+ylab("nu")
ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_weighted.png"), width=10, height=7.5)

#nu: correlations
nunowt = foreach (i=subs,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
    load(paste0("data/caulo_NcoI_150k_",i,"pcnowt_lambda",j,"_csnorm_optimized.RData"))
    data.table(method="as-is ", sub=i, lambda=j, logp=cs@par$logp, log_nu=cs@par$log_nu)
  }
nu = foreach (i=subs,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
  load(paste0("data/caulo_NcoI_150k_",i,"pc_lambda",j,"_csnorm_optimized.RData"))
  data.table(method="weighted ", sub=i, lambda=j, logp=cs@par$logp, log_nu=cs@par$log_nu)
}
nunowt[,log_nu_ref:=nu[sub==100,log_nu],by=sub]
nu[,log_nu_ref:=nu[sub==100,log_nu],by=sub]
#
ggplot(rbind(nu,nunowt)[,.(correlation=cor(log_nu,log_nu_ref)),by=c("method","sub")])+
  geom_point(aes(sub,correlation^2,colour=method))+geom_line(aes(sub,correlation^2,colour=method))+
  xlab("percent of counts used for fitting")+ylab("R-squared")
ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_correlation.png"), width=10, height=7.5)


#nu: local analyses
nu = foreach (i=subs,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
    load(paste0("data/caulo_NcoI_150k_",i,"pc_lambda",j,"_csnorm_optimized.RData"))
    #number of counts per cut site
    sums=rbind(cs@par$counts.sub[,.N,by=id1][,.(id=id1,N)],cs@par$counts.sub[,.N,by=id2][,.(id=id2,N)])[,.(N=sum(N)),keyby=id]
    sums=merge(sums,cs@biases[,.(id)],by="id",all.y=T)
    sums[is.na(N),N:=0]
    data.table(method="weighted ", sub=i, lambda=j, logp=cs@par$logp, log_nu=cs@par$log_nu, ncounts=sums[,N])
}
nu[,log_nu_ref:=nu[sub==100,log_nu],by=sub]
nu[,Sres:=(log_nu-log_nu_ref)^2]
nu[,Stot:=log_nu^2]
nu[,Rsq:=1-sum(Sres)/sum(Stot),by=sub]
ggplot(nu)+geom_point(aes(ncounts,Sres/sum(Sres)*100))+ylim(0,2)+scale_x_continuous(breaks=seq(0,80,by = 10))+
  xlab("Number of counts per cut site used for fitting (all=80)") + ylab("Relative contribution to chi square (%)")
ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_chisq.png"), width=10, height=7.5)





### reproducibility: initial guess vs full estimation
nreads=c(5,10,20,30,40,50,60,75,90,100,125,150,175,200,300,400)
lambdas=c(0.01,0.1,1,10,100)

registerDoParallel(cores=30)
foreach (nread=nreads) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",nread,"k_lambda",lambda,"_csnorm_optimized.RData"))
  op=list(par=cs@par)
  #log_decay=cs@par$log_decay
  load(paste0("data/caulo_NcoI_150k_sub",nread,"k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  bf_per_kb=0.25
  dmin=1-0.01
  dmax=150000+0.01
  #initial guess
  #init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
  #  counts=cs@counts, biases=cs@biases, lambda=lambda,
  #  bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, verbose=T, iter=10000)))
  #op=list(par=init.op)
  for (i in 1:3) {
    a=system.time(output <- capture.output(op <- csnorm:::csnorm_simplified(
      model=csnorm:::stanmodels$simplified, biases=cs@biases, counts = cs@counts,
      dmin=dmin, dmax=dmax, log_decay=op$par$log_decay, log_nu=op$par$log_nu, log_delta=op$par$log_delta,
      groups=10, bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = T, init=op)))
  }
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  op$par$value=op$value
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  #cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  #cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  #cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_150k_sub",nread,"k_lambda",lambda,"_csnorm_optimized_init.RData"))
}



# generate plots
prefix="NcoI_150k"

outputs = foreach (i=nreads, .combine=rbind) %:% foreach (j=lambdas, .combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
  data.table(dset=paste("nreads",i,"lambda",j),size=i, lambda=j, logp=cs@par$value,
             runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, out=tail(cs@par$output,1))
}
outputs
summary(lm(data=outputs, log(runtime)~log(size)+log(lambda)))
ggplot(outputs)+geom_point(aes(size,runtime,colour=log(lambda)))+scale_x_log10()+scale_y_log10()+
  ylab("run time (s)")+xlab("dataset size")+ggtitle("t ~ size^-0.36")
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_runtime.png"), width=10, height=7.5)

#nu: initial with all lambdas
registerDoParallel(cores=30)
nuinit = foreach (i=nreads,.combine=rbind)  %:% foreach (j=lambdas,.combine=rbind) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                     bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                       ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                          dset=paste("init nreads",i,"lambda",j),size=i, lambda=j,cat="init",logp=cs@par$value)]
}
ggplot(nuinit)+geom_line(aes(pos,exp(log_nu)))+facet_grid(lambda~size, labeller=label_both)+ylim(0,2)+ylab("nu")

#nu: simplified with all lambdas
registerDoParallel(cores=30)
nusimplified = foreach (i=nreads,.combine=rbind)  %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                        dset=paste("simplified nreads",i,"lambda",j),size=i, lambda=j,cat="simplified",logp=cs@par$value)]
}
nusimplified[,rank:=factor(frank(-logp,ties.method="dense")),by=size]
ggplot(nusimplified)+geom_line(aes(pos,exp(log_nu),colour=rank))+facet_grid(lambda~size, labeller=label_both)+ylim(0,2)+ylab("nu")

#nu: full with all lambdas
registerDoParallel(cores=30)
nufull = foreach (i=nreads,.combine=rbind)  %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                        dset=paste("full nreads",i,"lambda",j),size=i, lambda=j,cat="full",logp=cs@par$value)]
}
nufull[,rank:=factor(frank(-logp,ties.method="dense")),by=size]
ggplot(nufull)+geom_line(aes(pos,exp(log_nu),colour=rank))+facet_grid(lambda~size, labeller=label_both)+ylim(0,2)+ylab("nu")



#nu: simplified and full compared, all lambdas
registerDoParallel(cores=30)
nu=rbind(foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                     bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                       ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                          dset=paste("simplified nreads",i,"lambda",j),size=i, lambda=j,cat="simplified",logp=cs@par$value)]
  },
  foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                     bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                       ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                          dset=paste("full nreads",i,"lambda",j),size=i, lambda=j,cat="full",logp=cs@par$value)]
  })
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=cat))+facet_grid(lambda~size, labeller=label_both)+ylim(0,2)+ylab("nu")
#ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_weighted.png"), width=10, height=7.5)


#nu: simplified and full compared, optimal lambdas
registerDoParallel(cores=30)
nusimplified = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                        dset=paste("simplified nreads",i,"lambda",j),size=i, lambda=j,cat="simplified",logp=cs@par$value)]
}
nufull = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                        dset=paste("full nreads",i,"lambda",j),size=i, lambda=j,cat="full",logp=cs@par$value)]
}
nu=rbind(nusimplified,nufull)
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=cat))+facet_wrap(~size, labeller=label_both)+ylim(0,2)+ylab("nu")
#ggsave(filename=paste0("images/caulo_NcoI_150k_sub_reproducibility_nu_weighted.png"), width=10, height=7.5)


#diagonal decay
registerDoParallel(cores=30)
dsimplified = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
               size=i,lambda=j,key="dist")
  }
dfull = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
               size=i,lambda=j,key="dist")
  }
decay=rbind(dsimplified,dfull)
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_wrap(~size, labeller=label_both)+
  ylab("decay")+scale_y_log10()+scale_x_log10()







