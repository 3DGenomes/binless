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

bf_per_kb=signif(10^(seq(-1,1,length.out = 9))*1/4,digits=2)


### normalize different datasets on a single CPU
registerDoParallel(cores=30)
foreach (bpk=bf_per_kb) %dopar% {
  load(paste0("data/caulo_NcoI_150k_csnorm.RData"))
  bf_per_decade=5
  dmin=1-0.01
  dmax=150000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
                                         counts=cs@counts, biases=cs@biases,
                                         bf_per_kb=bpk, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    model=csnorm:::stanmodels$fit, biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
    bf_per_kb=bpk, bf_per_decade=bf_per_decade, iter=100000, verbose = F, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bpk, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_150k_bfpkb",bpk,"_csnorm_optimized.RData"))
}


### generate plots
prefix="NcoI_150k"

#outputs and runtime
outputs = foreach (i=bf_per_kb,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_csnorm_optimized.RData"))
  data.table(dset=i,out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime, iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta)
}
outputs
summary(lm(data=outputs, log(runtime)~log(dset)))
ggplot(outputs)+geom_point(aes(dset,runtime))+scale_x_log10()+scale_y_log10()+ylab("run time (s)")+xlab("basis functions per kb")+
  ggtitle("t ~ bf^0.56")
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_runtime.png"), width=10, height=7.5)


#nu and delta: init
nu = foreach (i=bf_per_kb,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=i)]
}
nu[,dset:=ordered(dset,levels=bf_per_kb)]
ggplot(nu[,.(dset,pos,nu=exp(log_nu))])+ylim(0,2)+
  geom_line(aes(pos,nu),colour="black")+facet_wrap(~dset)

#nu and delta: plots and correlation
nu = foreach (i=bf_per_kb,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_150k_bfpkb",i,"_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=i)]
}
#nu
ggplot(nu[,.(dset,pos,nu=exp(log_nu))])+geom_line(aes(pos,nu),colour="black")+facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_nu.png"), width=10, height=7.5)
#nu*delta
ggplot(nu[,.(dset,pos,delta=exp(log_nu+log_delta))])+
  geom_line(aes(pos,delta),colour="black")+facet_wrap(~dset)+ylim(0,2)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_delta.png"), width=10, height=7.5)
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


