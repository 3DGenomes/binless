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



### reproducibility: initial guess vs full estimation
nreads=c(5,10,20,30,40,50,60,75,90,100,125,150,175,200,300,400)
lambdas=c(0.01,0.1,1,10,100)

registerDoParallel(cores=30)
foreach (nread=nreads) %:% foreach (lambda=lambdas) %dopar% {
  #load(paste0("data/caulo_NcoI_150k_sub",nread,"k_lambda",lambda,"_csnorm_optimized.RData"))
  #init.op=list(par=cs@par)
  #log_decay=cs@par$log_decay
  load(paste0("data/caulo_NcoI_150k_sub",nread,"k_csnorm.RData"))
  #dmin=1-0.01
  #dmax=150000+0.01
  cs = run_gibbs(cs, design=NULL, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambda=lambda,
                       ngibbs = 1, iter=100000, fit.decay=T, fit.genomic=F)
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
#ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_runtime.png"), width=10, height=7.5)

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
ggsave(filename="gibbs_full2_nu.png", width=10, height=7.5)

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
ggsave(filename="gibbs_full2_nuopt.png", width=10, height=7.5)



#diagonal decay: simplified and full compared, all lambdas
registerDoParallel(cores=30)
decay=rbind(foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
             size=i,lambda=j,key="dist")
},
foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
             size=i,lambda=j,key="dist")
})
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_grid(lambda~size, labeller=label_both)+scale_y_log10()+scale_x_log10()

#diagonal decay: optimal lambdas
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
#ggsave(filename="gibbs_full2_decay.png", width=10, height=7.5)



#parameters
registerDoParallel(cores=4)
psimplified = foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
    data.table(logp=cs@par$value,cat="simplified",size=i,lambda=j,
               eC=cs@par$eC, eRJ=cs@par$eRJ, eDE=cs@par$eDE, lambda_diag=cs@par$lambda_diag,
               lambda_delta=cs@par$lambda_delta, lambda_nu=cs@par$lambda_nu, alpha=cs@par$alpha)
  }
pfull = foreach (i=nreads,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    data.table(logp=cs@par$value,cat="full",size=i,lambda=j,
               eC=cs@par$eC, eRJ=cs@par$eRJ, eDE=cs@par$eDE, lambda_diag=cs@par$lambda_diag,
               lambda_delta=cs@par$lambda_delta, lambda_nu=cs@par$lambda_nu, alpha=cs@par$alpha)
  }
params=rbind(psimplified,pfull)
setkey(params, size,lambda,cat)
ggplot(params,aes(size,exp(eC),colour=cat))+geom_point()+ylim(0,5)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_jitter()+scale_x_log10()
ggplot(params)+geom_point(aes(factor(paste(size,lambda)),lambda_nu,colour=cat))+ylim(0,25)
ggplot(params)+geom_point(aes(factor(paste(size,lambda)),lambda_delta,colour=cat))#+ylim(0,25)
ggplot(dcast(params[,.(cat,size,lambda,exp(eC))],size+lambda~cat)[,.(size,lambda,ratio=simplified/full)][order(ratio)])+
  geom_point(aes(size,ratio,colour=log(lambda)))+ylim(0,2)

