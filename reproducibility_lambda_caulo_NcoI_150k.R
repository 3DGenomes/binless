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
    biases=cs@biases, counts = cs@counts, dmin=dmin, dmax=dmax,
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







