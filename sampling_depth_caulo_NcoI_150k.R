library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

### read caulobacter dataset and normalize a 150k square at different sampling depths

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

#produce subsampled datasets
foreach (nreads=c(5,10,20,30,40,50,60,75,90,100,125,150,175,200,300,400)) %do% {
  load("data/caulo_NcoI_150k_csdata_with_data.RData")
  data=cs@data[sample(.N,nreads*1000)]
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=4042929)
  cs = new("CSdata", info=cs@info, settings=cs@settings,
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  cs=merge_cs_norm_datasets(list(cs))
  save(cs, file=paste0("data/caulo_NcoI_150k_sub",nreads,"k_csnorm.RData"))
}

### normalize different datasets on a single CPU
registerDoParallel(cores=30)
foreach (nreads=c(5,10,20,30,40,50,60,75,90,100,125,150,175,200,300,400)) %:% foreach (lambda=c(0.01,0.1,1,10,100)) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",nreads,"k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  bf_per_kb=0.25
  dmin=1-0.01
  dmax=150000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
                                         counts=cs@counts, biases=cs@biases, lambda=lambda,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = F, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  op$par$value=op$value
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_150k_sub",nreads,"k_lambda",lambda,"_csnorm_optimized.RData"))
}


### generate plots
prefix="NcoI_150k"
nreads=c(5,10,20,30,40,50,60,75,90,100,125,150,175,200,300,400)
lambdas=c(0.01,0.1,1,10,100)

outputs = foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  data.table(dset=paste("bf",i,"lam",j),size=i,lambda=j, logp=cs@par$value,
             runtime=cs@par$runtime+cs@par$init$runtime, iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta)
}
outputs
summary(lm(data=outputs, log(runtime)~log(size)+log(lambda)))
ggplot(outputs)+geom_point(aes(size,runtime,colour=log(lambda)))+scale_x_log10()+scale_y_log10()+
  ylab("run time (s)")+xlab("dataset size")+ggtitle("t ~ size^-0.36")
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_runtime.png"), width=10, height=7.5)

#lFC
lFC = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=paste("bf",i,"lam",j),size=i,lambda=j, logp=cs@par$value, lFC)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
lFC=rbind(lFC,get_cs_binned(cs,1,"CS")[,.(dset="all",size="all", lambda=1, logp=-1, lFC)])
lFC[,dset:=ordered(size,levels=c("all",nreads))]
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %do% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    get_cs_binned(cs,1,"CS")[,.(dset=paste("bf",i,"lam",j),size=i,lambda=j, logp=cs@par$value,
                                begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
mat=rbind(mat,get_cs_binned(cs,1,"CS")[,.(dset="all",size="all",lambda=1, logp=-1,
                                          begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)])
mat[,dset:=ordered(size,levels=c("all",nreads))]
ggplot(mat[dset!=300])+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[dset!=300&is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_normalized.png"), width=10, height=7.5)

#nu and delta: init
registerDoParallel(cores=30)
nu.init = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,size=i,lambda=j, logp=cs@par$value)]
}
nuref=csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                       bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                         ,.(pos,log_nu_ref=log_nu,log_delta_ref=log_delta,
                                            dset="all",size="all",lambda=1, logp=-1)]
nu.init[,c("log_nu_ref","log_delta_ref"):=nuref[,.(log_nu_ref,log_delta_ref)],by=size]
nu.init[,dset:=ordered(size,levels=nreads)]
ggplot(nu.init[,.(dset,pos,nu=exp(log_nu),nuref=exp(log_nu_ref))])+
  geom_line(aes(pos,nu),colour="black")+geom_line(aes(pos,nuref),colour="red")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu_init.png"), width=10, height=7.5)
ggplot(nu.init[,.(dset,pos,delta=exp(log_delta),deltaref=exp(log_delta_ref))])+
  geom_line(aes(pos,delta),colour="black")+geom_line(aes(pos,deltaref),colour="red")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_delta_init.png"), width=10, height=7.5)


#nu and delta: plots and correlation
registerDoParallel(cores=30)
nu = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,size=i,lambda=j, logp=cs@par$value)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
nuref=csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                       bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                         ,.(pos,log_nu_ref=log_nu,log_delta_ref=log_delta,
                                            dset="all",size="all",lambda=1, logp=-1)]
nu[,c("log_nu_ref","log_delta_ref"):=nuref[,.(log_nu_ref,log_delta_ref)],by=size]
nu[,dset:=ordered(size,levels=nreads)]
ggplot(nu[,.(dset,pos,nu=exp(log_nu),nuref=exp(log_nu_ref))])+
  geom_line(aes(pos,nu),colour="black")+geom_line(aes(pos,nuref),colour="red")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu.png"), width=10, height=7.5)
#
ggplot(nu[,.(dset,pos,delta=exp(log_delta),deltaref=exp(log_delta_ref))])+
  geom_line(aes(pos,delta),colour="black")+geom_line(aes(pos,deltaref),colour="red")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_delta.png"), width=10, height=7.5)
#
ggplot(rbind(nu[,.(var="log_nu",correlation=cor(log_nu_ref,log_nu)),by=dset],
             nu[,.(var="log_delta",correlation=cor(log_delta_ref,log_delta)),by=dset]),aes(dset,correlation,colour=var))+
  geom_point()+ylim(-0.02,1)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu_delta_correlation.png"), width=10, height=7.5)
#
ggplot(rbind(nu[,.(var="log_nu",chisq=mean((log_nu_ref-log_nu)^2)),by=dset],
             nu[,.(var="log_delta",chisq=mean((log_delta_ref-log_delta)^2)),by=dset]),aes(dset,chisq,colour=var))+geom_point()
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu_delta_chisq.png"), width=10, height=7.5)
#
get_fpt=function(x,y){a=acf(y,plot=F,lag.max=10000); l=a$lag[a$acf<=0][1]; x[l+1]-x[1]}
fpts=nu[,.(fpt=get_fpt(pos,exp(log_nu+log_delta)),fptref=get_fpt(pos,exp(log_nu_ref+log_delta_ref))),by=dset]
ggplot(fpts)+geom_point(aes(dset,fpt/1000))+geom_hline(aes(yintercept=fptref/1000))+
  ylim(0,60)+ylab("first passage time of nu x delta (in kb)")
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu_delta_fpt.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=nreads,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value, dset=i,key="dist")
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
decay=rbind(decay,data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=-1, dset="all",key="dist"))
decay[,dset:=ordered(dset,levels=c("all",nreads))]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_decay.png"), width=10, height=7.5)


