library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)

### take less reads
subs = round(100/(2**(0:7)), digits = 2)
foreach (sub=subs) %do% {
  load("data/rao_HiCall_GM12878_SELP_150k_csdata_with_data.RData")
  cs=merge_cs_norm_datasets(list(csd), different.decays="none")
  cs=subsample_csnorm(cs, subsampling.pc = sub)
  save(cs, file=paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csnorm.RData"))
}


#run normalizations
registerDoParallel(cores=8)
runs = foreach (sub=subs) %dopar% {
  load(paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csnorm.RData"))
  cs = normalize_binless(cs)
  save(cs,file=paste0("data/rao_HiCall_GM12878_SELP_150k_",sub,"pc_csnorm_optimized.RData"))
  cs
}

#check convergence
sapply(runs,function(x){has_converged(x)})

#iota and rho
iota = foreach(sub=subs,cs=runs,.combine=rbind) %dopar% {
  points=data.table(type="pts",position=cs@biases[,pos],iota=exp(cs@par$log_iota),rho=exp(cs@par$log_rho))
  biases=generate_genomic_biases(cs)[,.(type="fun",position=pos,iota=exp(log_iota),rho=exp(log_rho))]
  ret=rbind(points,biases)
  ret[,method:=sub]
  ret
}

ggplot()+facet_grid(method~.)+scale_y_log10()+coord_cartesian(ylim=c(0.1,10))+
  geom_line(aes(position,iota),data=iota[type=="fun"])+#xlim(1000000,1100000)+
  geom_vline(aes(xintercept=position),data=iota[type=="pts"],alpha=0.1)
ggsave(filename = "images/rao_HiCall_SELP_150k_subs_iota_bias.pdf", width=10, height=7)


#decay
decay = foreach(j=subs,cs=runs,.combine=rbind) %dopar% {
  if ("decay" %in% names(cs@par$decay)) {
    cs@par$decay[,.(subsampling=j,distance=dist,decay)]
  } else {
    cs@par$decay[,.(subsampling=j,distance,decay=exp(kappa-cs@par$eC))]
  }
}
ggplot(decay[,.SD[sample(.N,min(.N,100000))],by=subsampling])+scale_colour_gradient(trans="log10")+
  geom_line(aes(distance,decay,colour=subsampling,group=subsampling))+scale_x_log10()+scale_y_log10()
ggsave(filename="images/rao_HiCall_SELP_150k_subs_diagonal_decay.pdf", width=10, height=7)

