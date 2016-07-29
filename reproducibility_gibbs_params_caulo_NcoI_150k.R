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
ngroups=c(3,5,7,10,25,50,79)
lambdas=c(0.01,0.1,1,10,100)

registerDoParallel(cores=30)
foreach (group=ngroups) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_csnorm.RData"))
  cs = run_gibbs(cs, design=NULL, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=group, lambda=lambda,
                       ngibbs = 2, iter=100000, fit.decay=T, fit.genomic=T)
  save(cs, file=paste0("data/caulo_NcoI_150k_sub400k_lambda",lambda,"_gibbs1_",group,"groups_csnorm_optimized.RData"))
}



# generate plots
prefix="NcoI_150k"

outputs = foreach (i=ngroups, .combine=rbind) %:% foreach (j=lambdas, .combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
  data.table(dset=paste("ngroups",i,"lambda",j),ngroups=i, lambda=j, logp=cs@par$value,
             runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, out=tail(cs@par$output,1))
}
outputs
summary(lm(data=outputs, log(runtime)~log(ngroups)+log(lambda)))
ggplot(outputs)+geom_point(aes(ngroups,runtime,colour=log(lambda)))+scale_x_log10()+scale_y_log10()+
  ylab("run time (s)")+xlab("dataset ngroups")+ggtitle("t ~ ngroups^-0.36")
#ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_runtime.png"), width=10, height=7.5)

#nu: simplified and full compared, all lambdas
registerDoParallel(cores=30)
nu=rbind(foreach (i=ngroups,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
  csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                     bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                       ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                          dset=paste("simplified ngroups",i,"lambda",j),ngroups=i, lambda=j,cat="simplified",logp=cs@par$value)]
  },
  foreach (i=ngroups,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_csnorm_optimized.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                     bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                       ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                          dset=paste("full ngroups",i,"lambda",j),ngroups=i, lambda=j,cat="full",logp=cs@par$value)]
  })
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=cat))+facet_grid(lambda~ngroups, labeller=label_both)+ylim(0,2)+ylab("nu")
#ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_nuall.png", width=10, height=7.5)

#nu: simplified and full compared, optimal lambdas
registerDoParallel(cores=30)
nusimplified = foreach (i=ngroups,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                        dset=paste("simplified ngroups",i,"lambda",j),ngroups=i, lambda=j,cat="simplified",logp=cs@par$value)]
}
nufull = foreach (i=ngroups,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_csnorm_optimized.RData"))
    csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,lambda_nu=cs@par$lambda_nu,
                                        dset=paste("full ngroups",i,"lambda",j),ngroups=i, lambda=j,cat="full",logp=cs@par$value)]
}
nu=rbind(nusimplified,nufull)
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=cat))+facet_wrap(~ngroups, labeller=label_both)+ylim(0,2)+ylab("nu")
#ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_nu.png", width=10, height=7.5)



#diagonal decay: simplified and full compared, all lambdas
registerDoParallel(cores=30)
decay=rbind(foreach (i=ngroups,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
             ngroups=i,lambda=j,key="dist")
},
foreach (i=ngroups,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_csnorm_optimized.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
             ngroups=i,lambda=j,key="dist")
})
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_grid(lambda~ngroups, labeller=label_both)+scale_y_log10()+scale_x_log10()
#ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_decayall.png", width=10, height=7.5)

#diagonal decay: optimal lambdas
registerDoParallel(cores=30)
dsimplified = foreach (i=ngroups,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
               ngroups=i,lambda=j,key="dist")
  }
dfull = foreach (i=ngroups,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
               ngroups=i,lambda=j,key="dist")
  }
decay=rbind(dsimplified,dfull)
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_wrap(~ngroups, labeller=label_both)+
  ylab("decay")+scale_y_log10()+scale_x_log10()
#ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_decay.png", width=10, height=7.5)



#parameters
registerDoParallel(cores=30)
pinit = foreach (i=ngroups,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
  data.table(logp=NA,cat="init",ngroups=i,lambda=j,
             eC=cs@par$init$eC, eRJ=cs@par$init$eRJ, eDE=cs@par$init$eDE, lambda_diag=cs@par$init$lambda_diag,
             lambda_delta=cs@par$init$lambda_delta, lambda_nu=cs@par$init$lambda_nu, alpha=cs@par$init$alpha)
}
psimplified = foreach (i=ngroups,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_gibbs1_",i,"groups_csnorm_optimized.RData"))
  data.table(logp=cs@par$value,cat="simplified",ngroups=i,lambda=j,
               eC=cs@par$eC, eRJ=cs@par$eRJ, eDE=cs@par$eDE, lambda_diag=cs@par$lambda_diag,
               lambda_delta=cs@par$lambda_delta, lambda_nu=cs@par$lambda_nu, alpha=cs@par$alpha)
  }
pfull =  {
  load(paste0("data/caulo_NcoI_150k_sub400k_lambda",j,"_csnorm_optimized.RData"))
  data.table(logp=cs@par$value,cat="full",ngroups=NA,lambda=j,
               eC=cs@par$eC, eRJ=cs@par$eRJ, eDE=cs@par$eDE, lambda_diag=cs@par$lambda_diag,
               lambda_delta=cs@par$lambda_delta, lambda_nu=cs@par$lambda_nu, alpha=cs@par$alpha)
  }
params=rbind(pinit,psimplified,pfull)
setkey(params, ngroups,lambda,cat)
params=params[cat!="init"]
ggplot(params,aes(ngroups,exp(eC),colour=cat))+geom_point()+ylim(0,5)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_jitter()+scale_x_log10()
ggplot(params)+geom_point(aes(factor(paste(ngroups,lambda)),lambda_nu,colour=cat))+ylim(0,25)
ggplot(params)+geom_point(aes(factor(paste(ngroups,lambda)),lambda_delta,colour=cat))#+ylim(0,25)
ggplot(dcast(params[,.(cat,ngroups,lambda,exp(eC))],ngroups+lambda~cat)[,.(ngroups,lambda,ratio=simplified/full)][order(ratio)])+
  geom_point(aes(ngroups,ratio,colour=log(lambda)))+ylim(0,2)
dcast(params[,.(cat,ngroups,lambda,log(lambda_nu))],ngroups+lambda~cat)



