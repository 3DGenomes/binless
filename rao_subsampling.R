library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")



### take less reads
foreach (sub=c(1,10,20,30,40,50,60,70,80,90,100)) %do% {
  load("data/rao_HiCall_GM12878_SELP_150k_csdata_with_data.RData")
  cs=merge_cs_norm_datasets(list(csd), different.decays="none")
  cs=subsample_csnorm(cs, subsampling.pc = sub)
  save(cs, file=paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csnorm.RData"))
}


#approximate run
foreach (sub=c(1,10,20,30,40,50,60,70,80,90,100)) %dopar% {
  load(paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csdata.RData"))
  cs=merge_cs_norm_datasets(list(csd), different.decays="none")
  cs = run_gauss(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=10, ngibbs = 40, iter=100000, init_alpha=1e-7, ncounts = 1000000)
  save(cs,file=paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csnorm_optimized.RData"))
}



#plots

dsets = foreach (sub=c(1,10,20,30,40,50,60,70,80,90,100),.combine=c) %do% 
  paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csnorm_optimized.RData")
names = c(1,10,20,30,40,50,60,70,80,90,100)

dsets = c("data/rao_HiCall_GM12878_SELP_150k_1pc_csnorm_optimized.RData",
          "data/rao_HiCall_GM12878_SELP_150k_10pc_csnorm_optimized.RData",
          "data/rao_HiCall_GM12878_SELP_150k_10pc_csnorm_optimized_init1pc.RData",
          "data/rao_HiCall_GM12878_SELP_150k_20pc_csnorm_optimized.RData",
          "data/rao_HiCall_GM12878_SELP_150k_20pc_csnorm_optimized_init1pc.RData")
names=c("1","10","10i","20","20i")

#diagnostic plots
registerDoParallel(cores=30)
foreach(i=dsets,j=names) %dopar% {
  load(i)
  a=plot_diagnostics(cs)
  ggsave(a[[1]], filename=paste0("diag_",j,"pc.pdf"), width=10, height=8)
  ggsave(a[[2]], filename=paste0("diag2_",j,"pc.pdf"), width=10, height=8)
}

#iota and rho
iota = foreach(i=dsets,j=names,.combine=rbind) %dopar% {
  load(i)
  points=data.table(type="pts",pos=cs@biases[,pos],iota=exp(cs@par$log_iota),rho=exp(cs@par$log_rho))
  biases=get_genomic_biases(cs)[,.(type="fun",pos,iota=exp(log_iota),rho=exp(log_rho))]
  ret=rbind(points,biases)
  ret[,method:=j]
  ret
}

ggplot()+facet_grid(method~.)+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,iota),data=iota[type=="fun"])+#xlim(1000000,1100000)+
  geom_vline(aes(xintercept=pos),data=iota[type=="pts"],alpha=0.1)
ggsave(filename = "images/rao_HiCall_SELP_150k_subs_iota_bias.pdf", width=10, height=7)


#decay
decay = foreach(i=dsets,j=names,.combine=rbind) %dopar% {
  load(i)
  if ("decay" %in% names(cs@par$decay)) {
    cs@par$decay[,.(subsampling=j,dist,decay)]
  } else {
    cs@par$decay[,.(subsampling=j,dist=distance,decay=exp(kappa-cs@par$eC))]
  }
}
ggplot(decay[,.SD[sample(.N,min(.N,100000))],by=subsampling])+
  geom_line(aes(dist,decay,colour=subsampling,group=subsampling))+scale_x_log10()+scale_y_log10()
ggsave(filename="images/rao_HiCall_SELP_150k_subs_diagonal_decay.pdf", width=10, height=7)

#parameters
registerDoParallel(cores=30)
params = foreach(i=dsets,j=names,.combine=rbind) %dopar% {
  load(i)
  value=get_exact_logp(cs)
  data.table(subsampling=j,`exp eRJ`=exp(cs@par$eRJ),`exp eDE`=exp(cs@par$eDE),`exp eC`=exp(cs@par$eC),
             alpha=cs@par$alpha,`log lambda_iota`=log(cs@par$lambda_iota), `log lambda_rho`=log(cs@par$lambda_rho),
             `log lambda_diag`=log(cs@par$lambda_diag),loglik=cs@par$value, logp=value, ngibbs=cs@diagnostics$params[,.N/3])
}
params
params[,ngibbs:=NULL]
ggplot(melt(params,id.vars = "subsampling"))+
  geom_line(aes(subsampling,value))+geom_point(aes(subsampling,value))+facet_wrap(~ variable, scales="free_y")
ggsave(filename="images/rao_HiCall_SELP_150k_subs_parameters.pdf", width=10, height=7)
