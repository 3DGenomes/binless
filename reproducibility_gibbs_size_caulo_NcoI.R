library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(Hmisc)

setwd("/home/yannick/simulations/cs_norm")

#zoom on a portion of the dataset
load("data/caulo_NcoI_all_csdata_with_data.RData")
cs.ori=cs

registerDoParallel(cores=30)
foreach (i=seq(1,4,by=0.5)) %dopar% {
  data=cs.ori@data[re.closest1<=i*1000000&re.closest2<=i*1000000]
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=4042929)
  csd = new("CSdata", info=cs@info, settings=cs@settings,
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  cs = merge_cs_norm_datasets(list(csd))
  save(cs, file=paste0("data/caulo_NcoI_",i,"M_csnorm.RData"))
}


### reproducibility: initial guess vs full estimation
sizes=seq(1,4,by=0.5)
lambdas=c(0.01,0.05,0.1,0.5,1)

registerDoParallel(cores=30)
foreach (size=sizes) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_NcoI_",size,"M_csnorm.RData"))
  cs = run_gibbs(cs, design=NULL, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambda=lambda,
                       ngibbs = 1, iter=100000)
  save(cs, file=paste0("data/caulo_NcoI_",size,"M_lambda",lambda,"_gibbs1_csnorm_optimized.RData"))
}



# generate plots
outputs = foreach (i=sizes, .combine=rbind) %:% foreach (j=lambdas, .combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"M_lambda",j,"_gibbs1_csnorm_optimized.RData"))
  data.table(dset=paste("sizes",i,"lambda",j),size=i, lambda=j, logp=cs@par$value,
             runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, out=tail(cs@par$output,1))
}
outputs
summary(lm(data=outputs, log(runtime)~log(size)+log(lambda)))
ggplot(outputs)+geom_point(aes(size,runtime,colour=log(lambda)))+scale_x_log10()+scale_y_log10()+
  ylab("run time (s)")+xlab("dataset size")+ggtitle("t ~ size^-0.36")

#nu and delta: plots and correlation
registerDoParallel(cores=30)
nu = foreach (i=sizes, .combine=rbind) %:% foreach (j=lambdas, .combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"M_lambda",j,"_gibbs1_csnorm_optimized.RData"))
  rbind(csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,size=i, lambda=j, status="init")],
        csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                         bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                           ,.(pos,log_nu,log_delta,size=i, lambda=j, status=strsplit(tail(cs@par$output,1)," ")[[1]][3])])
}
ggplot(nu)+scale_y_log10(limits=c(0.1,10))+geom_line(aes(pos,exp(log_nu),colour=status))+facet_grid(size~lambda)+ylab("nu")
ggplot(nu)+scale_y_log10(limits=c(0.1,10))+geom_line(aes(pos,exp(log_delta),colour=status))+facet_grid(size~lambda)+ylab("delta")
ggplot(nu[pos>150000&pos<300000])+scale_y_log10()+geom_line(aes(pos,exp(log_nu),colour=status))+facet_grid(size~lambda)+ylab("nu")



#diagonal decay: simplified and full compared, all lambdas
registerDoParallel(cores=30)
decay=rbind(foreach (i=sizes,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
             size=i,lambda=j,key="dist")
},
foreach (i=sizes,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
             size=i,lambda=j,key="dist")
})
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_grid(lambda~size, labeller=label_both)+scale_y_log10()+scale_x_log10()
ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_decayall.png", width=10, height=7.5)

#diagonal decay: optimal lambdas
registerDoParallel(cores=30)
dsimplified = foreach (i=sizes,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
               size=i,lambda=j,key="dist")
  }
dfull = foreach (i=sizes,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
               size=i,lambda=j,key="dist")
  }
decay=rbind(dsimplified,dfull)
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_wrap(~size, labeller=label_both)+
  ylab("decay")+scale_y_log10()+scale_x_log10()
ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_decay.png", width=10, height=7.5)



#parameters
registerDoParallel(cores=30)
params = foreach (i=sizes,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"M_lambda",j,"_gibbs1_csnorm_optimized.RData"))
  rbind(data.table(logp=NA,cat="init",size=i,lambda=j,status=NA,
             eC=cs@par$init$eC, eRJ=cs@par$init$eRJ, eDE=cs@par$init$eDE, lambda_diag=cs@par$init$lambda_diag,
             lambda_delta=cs@par$init$lambda_delta, lambda_nu=cs@par$init$lambda_nu, alpha=cs@par$init$alpha),
        data.table(logp=cs@par$value,cat="gibbs",size=i,lambda=j,status=strsplit(tail(cs@par$output,1)," ")[[1]][3],
               eC=cs@par$eC, eRJ=cs@par$eRJ, eDE=cs@par$eDE, lambda_diag=cs@par$lambda_diag,
               lambda_delta=cs@par$lambda_delta, lambda_nu=cs@par$lambda_nu, alpha=cs@par$alpha))
}
setkey(params, size,lambda,cat)
ggplot(params,aes(size,eC,colour=cat))+geom_point()+ylim(-5,20)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_jitter()#+scale_x_log10()
ggplot(params)+geom_point(aes(factor(paste(size,lambda)),lambda_nu,colour=cat))+scale_y_log10()
ggplot(params)+geom_point(aes(factor(paste(size,lambda)),lambda_delta,colour=cat))#+ylim(0,25)
ggplot(dcast(params[,.(cat,size,lambda,eC)],size+lambda~cat)[,.(size,lambda,ratio=init/gibbs)][order(ratio)])+
  geom_point(aes(size,ratio,colour=log(lambda)))+ylim(0,2)
dcast(params[,.(cat,size,lambda,log(lambda_nu))],size+lambda~cat)




load(paste0("data/caulo_NcoI_3M_csnorm.RData"))
cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases)
cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
dmin=0.99
dmax=cs@settings$circularize/2+0.01
cs@settings$dmin=dmin
cs@settings$dmax=dmax
lambda=0.1
#init.a=system.time(init.output <- capture.output(init.par <- run_split_parallel_initial_guess(
#  counts=cs@counts, biases=cs@biases, lambda=0.1, verbose=T,
#  bf_per_kb=0.25, dmin=dmin, dmax=dmax, bf_per_decade=5, iter=1000000)))
init.a=system.time(init.output <- capture.output(init.op <- csnorm:::csnorm_simplified_guess(
  biases = cs@biases, counts = cs@counts, dmin = dmin, dmax = dmax, lambda=0.1,
  groups = 10, bf_per_kb = 0.25, bf_per_decade=5, iter = 1000000)))
init.par=init.op$par
a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_simplified_genomic(
  biases = cs@biases, counts = cs@counts,
  log_decay = init.par$log_decay, log_nu = init.par$log_nu, log_delta = init.par$log_delta,
  groups = 10, bf_per_kb = 0.25, iter = 100, init=init.par)))
cs@par=op.gen$par
nu=rbind(csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                       bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                         ,.(pos,log_nu,log_delta,size=i, lambda=j, status="init")],
      csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                       bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                         ,.(pos,log_nu,log_delta,size=i, lambda=j, status=strsplit(tail(cs@par$output,1)," ")[[1]][3])])

