library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

if (F) {
  csd=generate_fake_dataset(signal=F,eC=-4.4,replicate="1",condition="WT")
  save(csd,file="data/fake_replicate1_csdata.RData")
  biases.ref=csd@biases
  csd=generate_fake_dataset(signal=F,biases.ref=biases.ref,eC=-5,replicate="2",condition="WT")
  save(csd,file="data/fake_replicate2_csdata.RData")
  csd=generate_fake_dataset(signal=F,eC=-4,biases.ref=biases.ref,replicate="3",condition="WT")
  save(csd,file="data/fake_replicate3_csdata.RData")
  csd=generate_fake_dataset(signal=F,eC=-3.8,biases.ref=biases.ref,replicate="4",condition="WT")
  save(csd,file="data/fake_replicate4_csdata.RData")
  csd=generate_fake_dataset(signal=F,eC=-4.1,biases.ref=biases.ref,replicate="5",condition="WT")
  save(csd,file="data/fake_replicate5_csdata.RData")
  csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.3,replicate="1",condition="KO")
  save(csd,file="data/fake_signal_replicate1_csdata.RData")
  csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-6.1,replicate="2",condition="KO")
  save(csd,file="data/fake_signal_replicate2_csdata.RData")
  csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-3.7,replicate="3",condition="KO")
  save(csd,file="data/fake_signal_replicate3_csdata.RData")
  csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.2,replicate="4",condition="KO")
  save(csd,file="data/fake_signal_replicate4_csdata.RData")
  csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-5.2,replicate="5",condition="KO")
  save(csd,file="data/fake_signal_replicate5_csdata.RData")

  counts=csd@counts
  counts[,bin1:=round(pos1/10000)]
  counts[,bin2:=round(pos2/10000)]
  binned=counts[,.(count=sum(contact.far+contact.close+contact.up+contact.down)),by=c("bin1","bin2")]
  ggplot(binned)+geom_raster(aes(bin1,bin2,fill=log(count)))+scale_fill_gradient2()

}

args=commandArgs(trailingOnly=TRUE)
nbg=as.integer(args[1])
nsig=as.integer(args[2])
fit.signal=as.logical(args[3])
ncores=10

datasets=list()
if (nbg>=1) {
  for (i in 1:nbg) {
    load(paste0("data/fake_replicate",i,"_csdata.RData"))
    datasets=c(datasets,list(csd))
  }
}
if (nsig>=1) {
  for (i in 1:nsig) {
    load(paste0("data/fake_signal_replicate",i,"_csdata.RData"))
    datasets=c(datasets,list(csd))
  }
}

cs=merge_cs_norm_datasets(datasets, different.decays="none", dfuse=5,qmin=0.)

base.res=10000

cs = run_gauss(cs, restart=F, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10,
               ngibbs = 25, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
               ncounts = 1000000, ncores=5, base.res=base.res, fit.signal=fit.signal, fit.disp=T, fit.decay=T, fit.genomic=T)

save(cs,file=paste0("data/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_csnorm_optimized.RData"))

for (resolution in c(base.res,5000,20000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  if (nbg+nsig>1) {
    cs=detect_binned_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
    cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
  }
  save(cs,file=paste0("data/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_csnorm_optimized.RData"))
}


if (F) {
  
  load(paste0("data/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_csnorm_optimized.RData"))
  
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
  
  
  #observed
  resolution=10000
  mat=get_matrices(cs, resolution=resolution, group="all")
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=observed))+
    geom_raster(aes(begin2,begin1,fill=observed))+coord_fixed()+
    scale_fill_gradient(low="white",high="black",na.value="white",trans="log")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + #guides(colour=F) +
    theme_minimal()+ theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/fake_",nbg,"bg_",nsig,"sig_observed.png"),width=7,height=5)
  
  #normalized
  mat=get_matrices(cs, resolution=resolution, group="all")
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=normalized))+
    geom_raster(aes(begin2,begin1,fill=normalized))+coord_fixed()+facet_wrap(~name)+
    scale_fill_gradient(low="white",high="black",na.value="white",trans="log")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_void()+ theme(axis.title=element_blank(),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_binned_normalized.png"),width=7,height=5)
  
  load("data/fake_0bg_1sig_TRUE_csnorm_optimized.RData")
  csT=cs
  load("data/fake_0bg_1sig_FALSE_csnorm_optimized.RData")
  csF=cs
  matF=get_matrices(csF, resolution=resolution, group="all")
  matT=get_matrices(csT, resolution=resolution, group="all")
  ice=merge(matF[,.(name,begin1,begin2,observed)],iterative_normalization(matF,niterations = 10),by=c("name","begin1","begin2"))
  ice[,ice.10:=ice.10/70.5]
  mat=rbind(matF[,.(model="no signal",begin1,begin2,observed,normalized)],
            matT[,.(model="signal",begin1,begin2,observed,normalized)],
            ice[,.(model="ICE 10",begin1,begin2,observed,normalized=ice.10)])
  
  #normalized side by side
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=normalized))+
    geom_raster(aes(begin2,begin1,fill=normalized))+coord_fixed()+facet_wrap(~model)+
    scale_fill_gradient(low="white",high="black",na.value="white",trans="log")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_minimal()+ theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/fake_",nbg,"bg_",nsig,"sig_sidebyside_binned_normalized.png"),width=7,height=5)
  
  #column sums side by side
  colsums=rbind(mat[,.(model,begin1,begin2,normalized)],mat[begin1!=begin2,.(model,begin1=begin2,begin2=begin1,normalized)])[
    ,.(coverage=sum(normalized)),by=c("model","begin1")]
  ggplot(colsums)+geom_line(aes(begin1,coverage))+facet_wrap(~model)+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/fake_",nbg,"bg_",nsig,"sig_sidebyside_binned_coverage.png"),width=7,height=2)
  
  
}

