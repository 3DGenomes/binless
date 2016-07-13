library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

### read rao dataset

a=examine_dataset("/scratch/rao/mapped/HICall_both_filled_map_chrX_73780165-74230165.tsv",
                  skip="SRR",nrows=1000000)
csd=read_and_prepare("/scratch/rao/mapped/HICall_both_filled_map_chrX_73780165-74230165.tsv",
                     "data/rao_HiCall_chrX_450k", "WT", "1", skip="SRR", circularize=-1, dangling.L=c(-1,0,3,4,5,8),
                     dangling.R=c(4,3,0,-1,-2,-5), maxlen=1000, save.data=T)
cs=merge_cs_norm_datasets(list(csd))
save(cs, file="data/rao_HiCall_chrX_450k_csnorm.RData")

#zoom on a 150kb portion of the dataset
load("data/rao_HiCall_chrX_450k_csdata_with_data.RData")
data=cs@data[re.closest1>=73780159&re.closest1<=73780159+150000&re.closest2>=73780159&re.closest2<=73780159+150000]
cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
dset_statistics(cs_data$biases,cs_data$counts)
csd = new("CSdata", info=cs@info, settings=cs@settings,
          data=data, biases=cs_data$biases, counts=cs_data$counts)
cs=merge_cs_norm_datasets(list(csd))
save(cs, file="data/rao_HiCall_chrX_150k_csnorm.RData")

bf_per_kb=signif(10^(seq(-0.7,0.7,length.out = 6)),digits=2)
lambdas=c(0.01,0.1,1,10,100)


### normalize different datasets on a single CPU
registerDoParallel(cores=30)
foreach (bpk=bf_per_kb) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/rao_HiCall_chrX_150k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=abs(pos2-pos1)]
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
  op$par$logp=op$value
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bpk, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  #cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  #cs=postprocess(cs, resolution=2000, ncores=10, verbose=F)
  #cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/rao_HiCall_chrX_150k_bfpkb",bpk,"_lambda",lambda,"_csnorm_optimized.RData"))
}




### normalize different datasets on a single CPU
registerDoParallel(cores=30)
foreach (bpk=bf_per_kb) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/rao_HiCall_chrX_150k_bfpkb",bpk,"_lambda",lambda,"_csnorm_optimized.RData"))
  #cs@pred=csnorm_predict_all(cs,ncores=10,verbose=F)
  cs@binned=list()
  cs=postprocess(cs, resolution=5000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/rao_HiCall_chrX_150k_bfpkb",bpk,"_lambda",lambda,"_csnorm_optimized.RData"))
}




### generate plots
prefix="rao_HiCall_chrX_150k"

#outputs and runtime
outputs = foreach (i=bf_per_kb,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/rao_HiCall_chrX_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
  data.table(dset=paste("bf",i,"lam",j),bfpkb=i,lambda=j,
             out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime, iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$logp)
}
outputs
summary(lm(data=outputs, log(runtime)~log(bfpkb)+log(lambda)))
ggplot(outputs)+geom_point(aes(bfpkb,runtime,colour=log(lambda)))+
  ylab("run time (s)")+xlab("basis functions per kb")+ggtitle("t ~ bf^0.46")
ggsave(filename=paste0("images/",prefix,"_reproducibility_runtime.png"), width=10, height=7.5)


#nu and delta: init
registerDoParallel(cores=30)
get_fpt=function(x,y){a=acf(y,plot=F,lag.max=10000); l=a$lag[a$acf<=0][1]; x[l+1]-x[1]}
nu.init = foreach (i=bf_per_kb,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/rao_HiCall_chrX_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("bf",i,"lam",j),
                                        bfpkb=i,lambda=j,logp=-1,lambda_nu=cs@par$lambda_nu,
                                        fpt=get_fpt(pos,log_nu))]
}
nu.init[,rank:=frank(-logp,ties.method="dense"),by=bfpkb]
ggplot(nu.init)+geom_line(aes(pos,exp(log_nu)))+facet_grid(lambda~bfpkb)+ylim(0,2)+ scale_colour_gradientn(colours = rainbow(3))
ggsave(filename=paste0("images/",prefix,"_reproducibility_nu_init.png"), width=10, height=7.5)

#nu and delta: plots and correlation
registerDoParallel(cores=30)
nu = foreach (i=bf_per_kb,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/rao_HiCall_chrX_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=paste("bf",i,"lam",j),
                                        bfpkb=i,lambda=j,logp=cs@par$logp,lambda_nu=cs@par$lambda_nu)]
}
nu[,rank:=frank(-logp,ties.method="dense"),by=bfpkb]
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=as.factor(rank)))+facet_grid(lambda~bfpkb)+ylim(0,2)
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=log(lambda_nu)))+facet_grid(rank~bfpkb)+ylim(0,2)
ggplot(nu[rank==1])+geom_line(aes(pos,exp(log_nu),colour=log(lambda_nu)))+facet_wrap(~bfpkb)+ylim(0,2)
ggsave(filename=paste0("images/",prefix,"_reproducibility_nu.png"), width=10, height=7.5)
#nu*delta
ggplot(nu[,.(dset,pos,delta=exp(log_nu+log_delta))])+
  geom_line(aes(pos,delta),colour="black")+facet_wrap(~dset)+ylim(0,2)
#delta
ggplot(nu[,.(dset,pos,delta=exp(log_delta))])+
  geom_line(aes(pos,delta),colour="black")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/",prefix,"_reproducibility_delta.png"), width=10, height=7.5)
#fpt
get_fpt=function(x,y){a=acf(y,plot=F,lag.max=10000); l=a$lag[a$acf<=0][1]; x[l+1]-x[1]}
fpts=nu[,.(fpt=get_fpt(pos,exp(log_nu+log_delta))),by=dset]
ggplot(fpts)+geom_point(aes(4/dset,fpt/1000))+
  ylab("first passage time of nu x delta (in kb)")+xlab("span of one basis function (in kb)")+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/",prefix,"_reproducibility_nu_delta_fpt.png"), width=10, height=7.5)

#lFC
lFC = foreach (i=bf_per_kb,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
    load(paste0("data/rao_HiCall_chrX_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
    get_cs_binned(cs,1,"CS")[,.(dset=i,lFC,logp=cs@par$logp)]
}
lFC[,dset:=ordered(dset,levels=bf_per_kb)]
ggplot(lFC)+geom_violin(aes(dset,lFC))
ggsave(filename=paste0("images/",prefix,"_reproducibility_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=bf_per_kb,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
    load(paste0("data/rao_HiCall_chrX_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
    get_cs_binned(cs,1,"CS")[,.(dset=i,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected,
                                bfpkb=i,lambda=j,logp=cs@par$logp)]
  }
mat[,dset:=ordered(bfpkb,levels=bf_per_kb)]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/",prefix,"_reproducibility_normalized.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=bf_per_kb,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
    load(paste0("data/rao_HiCall_chrX_150k_bfpkb",i,"_lambda",j,"_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay), dset=i,key="dist", logp=cs@par$logp)
}
decay[,dset:=ordered(dset,levels=bf_per_kb)]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/",prefix,"_reproducibility_decay.png"), width=10, height=7.5)


