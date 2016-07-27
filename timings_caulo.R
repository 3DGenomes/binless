library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

### read caulobacter dataset and normalize a 150k square at different sampling depths

#zoom on a portion of the dataset
load("data/caulo_NcoI_all_csdata_with_data.RData")
cs.ori=cs

registerDoParallel(cores=30)
foreach (i=seq(150,1000,by=50)) %dopar% {
  data=cs.ori@data[re.closest1>=2000000&re.closest1<=2000000+i*1000&re.closest2>=2000000&re.closest2<=2000000+i*1000]
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=4042929)
  csd = new("CSdata", info=cs@info, settings=cs@settings,
           data=data, biases=cs_data$biases, counts=cs_data$counts)
  cs = merge_cs_norm_datasets(list(csd))
  save(cs, file=paste0("data/caulo_NcoI_",i,"k_csnorm.RData"))
}


### normalize different datasets on a single CPU
registerDoParallel(cores=30)
foreach (i=seq(150,500,by=50)) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  bf_per_kb=0.25
  dmin=1-0.01
  dmax=i*1000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
                                         counts=cs@counts, biases=cs@biases,
                                         bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade,
                                         verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = T, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  op$par$mem=as.integer(object.size(list(cs@biases,cs@counts)))
  op$par$nsteps=as.integer(strsplit(tail(output,3)[[1]],'\\s+')[[1]][2])
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_",i,"k_csnorm_optimized_bfpkb",bf_per_kb,".RData"))
}


### generate plots
sizes=seq(150,500,by=50)
bases=rep(bf_per_kb,length(sizes))
fnames=paste0("data/caulo_NcoI_",sizes,"k_csnorm_optimized_bfpkb",bases,".RData")
dnames=paste0(sizes,"k")
prefix="NcoI_bfpkb1"

#outputs and runtime
outputs = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  data.table(dset=d,size=s, bfpkb=b, out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, mem=cs@par$mem)
}
outputs
summary(lm(data=outputs, log(runtime)~log(size)))
ggplot(outputs)+geom_point(aes(sizes,runtime))+scale_x_log10()+scale_y_log10()+ylab("run time (s)")+xlab("genome size")+
  ggtitle("t ~ sz^2")
ggsave(filename=paste0("images/caulo_",prefix,"_timings_runtime.png"), width=10, height=7.5)


#nu and delta: init
nu.init = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=d,size=s, bfpkb=b)]
}
nu.init[,dset:=ordered(dset,levels=dnames)]
setkey(nu.init,dset,pos)
ggplot(nu.init[,.(dset,pos,nu=exp(log_nu))])+scale_y_log10(breaks=seq)+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
#
ggplot(nu.init[,.(dset,pos,nu=exp(log_nu))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_nu_init.png"), width=10, height=7.5)

#nu and delta: plots and correlation
nu = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=d,size=s, bfpkb=b)]
}
nu[,dset:=ordered(dset,levels=dnames)]
setkey(nu,dset,pos)
#nu
ggplot(nu[,.(dset,pos,nu=exp(log_nu))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_nu.png"), width=10, height=7.5)
#nu*delta
ggplot(nu[,.(dset,pos,nu=exp(log_nu+log_delta))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
#delta
ggplot(nu[,.(dset,pos,delta=exp(log_delta))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,delta),colour="black")+facet_grid(~dset)+coord_flip()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_delta.png"), width=10, height=7.5)
#fpt
get_fpt=function(x,y){a=acf(y,plot=F,lag.max=10000); l=a$lag[a$acf<=0][1]; x[l+1]-x[1]}
fpts=nu[,.(fpt=get_fpt(pos,exp(log_nu+log_delta))),by=size]
ggplot(fpts)+geom_point(aes(size,fpt/1000))+
  ylab("first passage time of nu x delta (in kb)")+xlab("genome size (in kb)")#+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_nu_delta_fpt.png"), width=10, height=7.5)

#lFC
lFC = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(lFC,dset=d,size=s, bfpkb=b)]
}
lFC[,dset:=ordered(dset,levels=dnames)]
ggplot(lFC)+geom_violin(aes(dset,lFC))
ggsave(filename=paste0("images/caulo_",prefix,"_timings_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(begin1,begin2,normalized,is.interaction,prob.observed.gt.expected,dset=d,size=s, bfpkb=b)]
}
mat[,dset:=ordered(dset,levels=dnames)]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_timings_normalized.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay), dset=d,size=s, bfpkb=b, key="dist")
}
decay[,dset:=ordered(dset,levels=dnames)]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_decay.png"), width=10, height=7.5)








### normalize different datasets on a single CPU
registerDoParallel(cores=30)
foreach (i=c(4,2,1,0.5,0.25,0.125,0.0625)) %dopar% {
  load(paste0("data/caulo_NcoI_350k_csnorm.RData"))
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  bf_per_decade=5
  bf_per_kb=i
  dmin=1-0.01
  dmax=350000+0.01
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade,
    verbose=F, iter=10000)))
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = T, init=init.op)))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  op$par$mem=as.integer(object.size(list(cs@biases,cs@counts)))
  op$par$nsteps=as.integer(strsplit(tail(output,3)[[1]],'\\s+')[[1]][2])
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=copy(csnorm_predict_all(cs,ncores=10,verbose=F))
  cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_350k_csnorm_optimized_bfpkb",bf_per_kb,".RData"))
}


### generate plots
bases=c(4,2,1,0.5,0.25,0.125,0.0625)
sizes=rep(350,length(bases))
fnames=paste0("data/caulo_NcoI_",sizes,"k_csnorm_optimized_bfpkb",bases,".RData")
dnames=paste0(bases,"bfpkb")
prefix="NcoI_bfpkb"

#outputs and runtime
outputs = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  data.table(dset=d,size=s, bfpkb=b, out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, mem=cs@par$mem)
}
outputs
summary(lm(data=outputs, log(runtime)~log(bfpkb)))
ggplot(outputs)+geom_point(aes(bfpkb,runtime))+scale_x_log10()+scale_y_log10()+ylab("run time (s)")+xlab("basis functions per kb")+
  ggtitle("t ~ bf")
ggsave(filename=paste0("images/caulo_",prefix,"_timings_runtime.png"), width=10, height=7.5)


#nu and delta: init
nu.init = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=d,size=s, bfpkb=b)]
}
nu.init[,dset:=ordered(dset,levels=dnames)]
setkey(nu.init,dset,pos)
ggplot(nu.init[,.(dset,pos,nu=exp(log_nu))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_nu_init.png"), width=10, height=7.5)

#nu and delta: plots and correlation
nu = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 100)[
                                     ,.(pos,log_nu,log_delta,dset=d,size=s, bfpkb=b)]
}
nu[,dset:=ordered(dset,levels=dnames)]
setkey(nu,dset,pos)
#nu
ggplot(nu[,.(dset,pos,nu=exp(log_nu))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_nu.png"), width=10, height=7.5)
#nu*delta
ggplot(nu[,.(dset,pos,nu=exp(log_nu+log_delta))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,nu),colour="black")+facet_grid(~dset)+coord_flip()
#delta
ggplot(nu[,.(dset,pos,delta=exp(log_delta))])+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,delta),colour="black")+facet_grid(~dset)+coord_flip()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_delta.png"), width=10, height=7.5)
#fpt
get_fpt=function(x,y){a=acf(y,plot=F,lag.max=10000); l=a$lag[a$acf<=0][1]; x[l+1]-x[1]}
fpts=nu[,.(fpt=get_fpt(pos,exp(log_nu+log_delta))),by=bfpkb]
ggplot(fpts)+geom_point(aes(bfpkb,fpt/1000))+
  ylab("first passage time of nu x delta (in kb)")+xlab("basis functions per kb")+scale_x_log10()#+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_nu_delta_fpt.png"), width=10, height=7.5)

#lFC
lFC = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(lFC,dset=d,size=s, bfpkb=b)]
}
lFC[,dset:=ordered(dset,levels=dnames)]
ggplot(lFC)+geom_violin(aes(dset,lFC))
ggsave(filename=paste0("images/caulo_",prefix,"_timings_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(begin1,begin2,normalized,is.interaction,prob.observed.gt.expected,dset=d,size=s, bfpkb=b)]
}
mat[,dset:=ordered(dset,levels=dnames)]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_timings_normalized.png"), width=10, height=7.5)

#diagonal decay
decay = foreach (i=fnames, d=dnames, s=sizes, b=bases, .combine=rbind) %do% {
  load(i)
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay), dset=d,size=s, bfpkb=b, key="dist")
}
decay[,dset:=ordered(dset,levels=dnames)]
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_timings_decay.png"), width=10, height=7.5)



