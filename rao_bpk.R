library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)
library(GGally)

setwd("/home/yannick/simulations/cs_norm")





#approximate run
foreach (sub=c(1,10,20,30,40,50,60,70,80,90,100)) %dopar% {
  load(paste0("data/rao_HiCall_GM12878_Peak1_450k_",sub,"pc_csdata.RData"))
  cs=merge_cs_norm_datasets(list(csd), different.decays="none")
  cs = normalize_binless(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=10, ngibbs = 40, iter=100000, init_alpha=1e-7, ncounts = 1000000)
  save(cs,file=paste0("data/rao_HiCall_GM12878_Peak1_450k_",sub,"pc_csnorm_optimized.RData"))
}



#plots

dsets = foreach (bpk=c(0.1,0.25,0.5,1,2,5,10),.combine=c) %do% 
  paste0("data/rao_HiCall_GM12878_Peak1_450k_bpk",bpk,"_csnorm_optimized.RData")
names = c(0.1,0.25,0.5,1,2,5,10)


#diagnostic plots
foreach(i=dsets,j=names) %do% {
  load(i)
  ggsave(cs@diagnostics$plot, filename=paste0("diag_bpk",j,".pdf"), width=10, height=8)
  ggsave(cs@diagnostics$plot2, filename=paste0("diag2_bpk",j,".pdf"), width=10, height=8)
}

#iota and rho
iota = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  points=data.table(type="pts",pos=cs@biases[,pos],iota=exp(cs@par$log_iota),rho=exp(cs@par$log_rho))
  biases=get_genomic_biases(cs)[,.(type="fun",pos,iota=exp(log_iota),rho=exp(log_rho))]
  ret=rbind(points,biases)
  ret[,method:=j]
  ret
}

ggplot()+facet_grid(method~.)+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,iota),data=iota[type=="fun"])
ggsave(filename = "images/rao_HiCall_Peak1_450k_bpk_iota_bias.pdf", width=10, height=7)
ggplot()+facet_grid(method~.)+scale_y_log10(limits=c(0.1,10))+
  geom_line(aes(pos,iota),data=iota[type=="pts"])

diota=dcast(iota[type=="pts",.(pos,iota,method)], pos~method, value.var="iota")
diota[,pos:=NULL]
pairs(diota,log="xy",xlim=c(0.1,10),ylim=c(0.1,10))


#decay
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  if ("decay" %in% names(cs@par$decay)) {
    cs@par$decay[,.(subsampling=j,dist,decay)]
  } else {
    cs@par$decay[,.(subsampling=j,dist=distance,decay=exp(kappa-cs@par$eC))]
  }
}
ggplot(decay[,.SD[sample(.N,min(.N,100000))],by=subsampling])+
  geom_line(aes(dist,decay,colour=subsampling,group=subsampling))+scale_x_log10()+scale_y_log10()
ggsave(filename="images/rao_HiCall_Peak1_450k_bpk_diagonal_decay.pdf", width=10, height=7)

#parameters
registerDoParallel(cores=30)
params = foreach(i=dsets,j=names,.combine=rbind) %dopar% {
  load(i)
  value=get_exact_logp(cs)
  data.table(subsampling=j,eRJ=cs@par$eRJ,eDE=cs@par$eDE,eC=cs@par$eC,alpha=cs@par$alpha,lambda_iota=cs@par$lambda_iota,
             lambda_rho=cs@par$lambda_rho,lambda_diag=cs@par$lambda_diag,logp=value)
}
params
ggplot(melt(params,id.vars = "subsampling"))+scale_x_log10()+
  geom_line(aes(subsampling,value))+geom_point(aes(subsampling,value))+facet_wrap(~ variable, scales="free_y")
ggsave(filename="images/rao_HiCall_Peak1_450k_bpk_parameters.pdf", width=10, height=7)
