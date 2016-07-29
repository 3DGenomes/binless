library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")



### Compare caulobacter normalized in parallel vs normalized in serial with subsampling


#normalize in parallel
load("data/caulo_NcoI_1000k_csnorm.RData")
cs@settings$circularize=-1
bf_per_kb=0.25
coverage=4
square.size=150000
cs=run_split_parallel(cs, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                      bf_per_decade=5, distance_bins_per_decade=100, lambdas=c(0.01,1,100),
                      verbose = F, iter=10000, ncores=30,
                      homogenize=F, outprefix="tmp/test")#, ops.count=ops.count, ops.bias=ops.bias)
cs=run_split_parallel_recovery(cs, "tmp/test", square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                               bf_per_decade=5, distance_bins_per_decade=100, lambdas=c(0.01,1,100), verbose = F,
                               iter=10000, ncores=30, homogenize=F)
cs@pred=csnorm_predict_all(cs, ncores=30)
cs=postprocess(cs, resolution=10000, ncores=30, verbose=F)
cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
save(cs, file="data/caulo_NcoI_1000k_csnorm_optimized.RData")


#normalize in serial
registerDoParallel(cores=30)
sub=10
foreach (lambda=c(0.01,1,100)) %dopar% {
  load("data/caulo_NcoI_1000k_csnorm.RData")
  bf_per_kb=0.25
  bf_per_decade=5
  dmin=1-0.01
  dmax=cs@settings$circularize/2+0.01
  cs@counts=fill_zeros(counts = cs@counts, biases = cs@biases)
  cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
  init.a=system.time(init.output <- capture.output(init.op <- csnorm:::run_split_parallel_initial_guess(
    counts=cs@counts, biases=cs@biases,
    bf_per_kb=bf_per_kb, dmin=dmin, dmax=dmax, bf_per_decade=bf_per_decade, lambda=lambda, verbose=F, iter=100000)))
  counts.sub=cs@counts[sample(.N,round(sub/100*.N))]
  a=system.time(output <- capture.output(op <- csnorm:::csnorm_fit(
    biases=cs@biases, counts = counts.sub, dmin=dmin, dmax=dmax,
    bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, iter=100000, verbose = T,
    init=init.op, weight=cs@counts[,.N]/counts.sub[,.N])))
  op$par$runtime=a[1]+a[4]
  op$par$output=output
  op$par$logp=op$value
  init.op$runtime=init.a[1]+init.a[4]
  init.op$output=init.output
  op$par$init=init.op
  op$par$counts.sub=counts.sub
  cs@par=op$par
  cs@settings = c(cs@settings, list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, dmin=dmin, dmax=dmax))
  cs@pred=csnorm_predict_all(cs,ncores=30,verbose=F)
  #cs=postprocess(cs, resolution=10000, ncores=10, verbose=F)
  #cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_1000k_",sub,"pc_lambda",lambda,"_csnorm_optimized.RData"))
}

registerDoParallel(cores=3)
foreach (lambda=c(0.01,1,100)) %dopar% {
  load(paste0("data/caulo_NcoI_1000k_",sub,"pc_lambda",lambda,"_csnorm_optimized.RData"))
  bf_per_kb=0.25
  bf_per_decade=5
  dmin=1-0.01
  dmax=cs@settings$circularize/2+0.01
  cs@settings$dmin=dmin
  cs@settings$dmax=dmax
  cs@pred=csnorm_predict_all(cs,ncores=30,verbose=F)
  cs=postprocess(cs, resolution=10000, ncores=30, verbose=F)
  cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
  save(cs, file=paste0("data/caulo_NcoI_1000k_",sub,"pc_lambda",lambda,"_csnorm_optimized.RData"))
}



#normalize with gibbs sampler
registerDoParallel(cores=3)
foreach (lambda=c(0.01,1,100)) %dopar% {
  load("data/caulo_NcoI_1000k_csnorm.RData")
  cs = run_gibbs(cs, design=NULL, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambda=lambda,
                 ngibbs = 1, iter=10000)
  save(cs, file=paste0("data/caulo_NcoI_1000k_gibbs1_lambda",lambda,"_csnorm_optimized.RData"))
}





### generate plots
prefix="NcoI_1000k"
lambdas=c(0.01,1,100)

#outputs and runtime: serial subsampled
sub=10
outputs = foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_1000k_10pc_lambda",j,"_csnorm_optimized.RData"))
  data.table(dset=j,
             out=tail(cs@par$output,1), runtime=cs@par$runtime+cs@par$init$runtime,
             iruntime=cs@par$init$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$logp)
}
outputs
#outputs and runtime: gibbs
outputs = foreach (j=lambdas,.combine=rbind) %do% {
  load(paste0("data/caulo_NcoI_1000k_gibbs1_lambda",j,"_csnorm_optimized.RData"))
  data.table(dset=j,
             out=tail(cs@par$output,1), runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, logp=cs@par$value)
}
outputs



fnames=c("data/caulo_NcoI_1000k_10pc_lambda100_csnorm_optimized.RData",
         "data/caulo_NcoI_1000k_gibbs1_lambda1_csnorm_optimized.RData",
         "data/caulo_NcoI_1000k_csnorm_optimized.RData")
dsets=c("serial 10%","gibbs", "parallel")


#lFC
lFC = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=j,lFC)]
}
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  get_cs_binned(cs,1,"CS")[,.(dset=j,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_normalized.png"), width=10, height=5)

#nu and delta correlation
nu = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],log_nu=cs@par$log_nu,log_delta=cs@par$log_delta,dset=j,key="pos")
}
ggplot(nu)+geom_line(aes(pos,log_nu,colour=dset))
ggplot(nu)+geom_line(aes(pos,log_nu,colour=dset))+xlim(2e6,2.15e6)
ggplot(nu)+geom_line(aes(pos,log_nu,colour=dset))+xlim(1e6,1.15e6)
nu[,pbin:=cut(pos,3)]
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=dset))+facet_wrap(~pbin,scales = "free_x", nrow=3)+
  scale_y_continuous(limits = c(0,2))+ylab("nu")
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_nu.png"), width=10, height=7.5)
ggplot(nu)+geom_line(aes(pos,exp(log_delta),colour=dset))+facet_wrap(~pbin,scales = "free_x", nrow=3)+
  scale_y_continuous(limits = c(0,2))+ylab("delta")

#diagonal decay
decay = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  cs@par$decay[,.(dist,decay,dset=j)]
}
decay[,decay:=decay/exp(mean(log(decay))),by=dset]
decay=decay[dist>100]
ggplot(decay)+geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/caulo_",prefix,"_reproducibility_decay.png"), width=10, height=7.5)
