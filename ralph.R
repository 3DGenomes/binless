library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)

setwd("/home/yannick/simulations/cs_norm")



if (FALSE) {
  a=examine_dataset("/scratch/ralph/HiC/B_rep1_Sox2.tsv", skip=0L,nrows=1000000, skip.fbm=F, read.len=75)
  
  foreach (cell=c("B","ES")) %:% foreach (replicate=c("1","2")) %:% foreach (locus=c("Sox2")) %do% {
    csd0=read_and_prepare(paste0("/scratch/ralph/HiC/3_Mapped/",cell,"_rep",replicate,"_",locus,".tsv"), locus=NULL,
                          paste0("data/ralph_",locus,"_",cell,"_rep",replicate), cell, replicate, enzyme="MboI",
                                name=paste0(locus,"_",cell,"_rep",replicate), circularize=-1,
                          dangling.L=c(0), dangling.R=c(3), maxlen=700, read.len=c(75,76), dmin=1000, save.data=T)
    NULL
  }
  
  load("data/ralph_Sox2_B_rep1_csdata_with_data.RData")
  plot_binned(csd@data, resolution=10000, b1=csd@data[,min(rbegin1)], e1=csd@data[,max(rbegin2)]) 
  plot_raw(csd@data, b1=csd@data[,min(rbegin1)+6500], e1=csd@data[,min(rbegin1)+8500]) 
  data=csd@data[rbegin1>min(rbegin1)+6500&rbegin2<min(rbegin1)+8500]
  data[,category:=grepl("[#~]",id)]
}




args=commandArgs(trailingOnly=TRUE)
sub=args[1]
bpk=30
dfuse=20 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
ncores=10


locus="Sox2"
cell="B"
load(paste0("data/ralph_",locus,"_",cell,"_rep1_csdata.RData"))
csd1=csd
cell="ES"
load(paste0("data/ralph_",locus,"_",cell,"_rep1_csdata.RData"))
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none", dfuse=dfuse)
cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 5, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
               ncounts = 1000000, ncores=ncores, base.res=10000, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T,
               signif.threshold=T)
save(cs,file=paste0("data/ralph_Sox2_B1_ES1_csnorm_optimized_base10k_cv.RData"))

cs = run_gauss(cs, restart=T, ngibbs = 15, ncores=ncores)
save(cs,file=paste0("data/ralph_Sox2_B1_ES1_csnorm_optimized_base10k_cv.RData"))

#load(paste0("data/ralph_Sox2_B1_ES1_csnorm_optimized_base10k_cv.RData"))

for (resolution in c(10000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores, signif.threshold=T)
  cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name], signif.threshold=T)
  save(cs,file=paste0("data/ralph_Sox2_B1_ES1_csnorm_optimized_base10k_cv.RData"))
}

if (F) {
  load("data/ralph_Sox2_csnorm_optimized_base10k.RData")
  
  csnorm:::has_converged(cs)
  cs@diagnostics$params[,sum(runtime)/3600]
  ggplot(cs@diagnostics$params[,.(step,leg,runtime)])+geom_line(aes(step,runtime,colour=leg))+scale_y_log10()
  
  plot_diagnostics(cs)$plot
  plot_diagnostics(cs)$plot2
  
  signals=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
    if ("signal" %in% cs@diagnostics$params[step==i,leg]) {
      sig=copy(cs@diagnostics$params[step==i&leg=="signal",signal][[1]])
      sig[,step:=i]
      sig
    }
  }
  
  ggplot(signals[name==name[1]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[name==name[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]-1])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_grid(step~ name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phihat))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
  
  #diagonal decay
  ggplot(cs@par$decay[name==name[1]&std<10])+geom_line(aes(distance,log_decay))+
    geom_line(aes(distance,log_decay+std),linetype="dotted")+geom_line(aes(distance,log_decay-std),linetype="dotted")+
    scale_x_log10()
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_decay.pdf",width=20,height=7)
  
  write.table(cs@par$decay[name==name[1]&std<10,.(distance,log_decay,std)],
              file="data/ralph_Sox2_B1_ES1_base10k_decay.txt", row.names=F, quote=F)
  
  
  resolution=10000
  mat.binned=get_matrices(cs, resolution=resolution, group="all")
  #observed
  ggplot(mat.binned)+geom_raster(aes(begin1,begin2,fill=log(observed)))+geom_raster(aes(begin2,begin1,fill=log(observed)))+facet_wrap(~name)+
    scale_fill_gradient(high="black", low="white", na.value = "white")+coord_fixed()
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_observed.pdf",width=20,height=7)
  
  mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi))+
    geom_raster(aes(begin2,begin1,fill=phi))+coord_fixed()+
    facet_wrap(~name)+scale_fill_gradient(high="black", low="white",na.value = "white")
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_binless.pdf",width=20,height=7)
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=beta))+
    geom_raster(aes(begin2,begin1,fill=beta))+coord_fixed()+
    facet_wrap(~name)+scale_fill_gradient(high="black", low="white",na.value = "white")
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_binless_nosignif.pdf",width=20,height=7)
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(weight)))+
    geom_raster(aes(begin2,begin1,fill=log(weight)))+coord_fixed()+
    facet_wrap(~name)+scale_fill_gradient(high="black", low="white",na.value = "white")
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_weight.pdf",width=20,height=7)
  
  mat.diff=get_interactions(cs, type="CSbdiff", resolution=resolution, group="all", ref=as.character(cs@experiments[1,name]))
  ggplot(mat.diff)+geom_raster(aes(begin1,begin2,fill=delta))+coord_fixed()+
    geom_raster(aes(begin2,begin1,fill=delta))+
    facet_wrap(~name)+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"),na.value="white")
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_bdiff.pdf",width=15,height=7)
  ggplot(mat.diff)+geom_raster(aes(begin1,begin2,fill=beta))+coord_fixed()+
    geom_raster(aes(begin2,begin1,fill=beta))+
    facet_wrap(~name)+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"),na.value="white")
  ggsave(filename="images/ralph_Sox2_B1_ES1_base10k_bdiff_nosignif.pdf",width=15,height=7)
  
  write.table(merge(merge(mat.binned[,.(name,begin1,begin2,observed)],
                          mat[,.(name,begin1,begin2,phi,beta.phi=beta,weight)]),
                    mat.diff[,.(name,begin1,begin2,delta,beta.delta=beta)]),
              file="data/ralph_Sox2_B1_ES1_base10k_normalized.txt", row.names=F, quote=F)
  
  
}  
  
