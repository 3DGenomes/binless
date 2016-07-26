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






