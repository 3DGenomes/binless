library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)

if (F) {
  csd=generate_fake_dataset(signal=F,eC=-4.4,replicate="1",condition="WT")
  save(csd,file="data/fake_replicate1_csdata.RData")
  biases.ref=csd@biases
  # csd=generate_fake_dataset(signal=F,biases.ref=biases.ref,eC=-5,replicate="2",condition="WT")
  # save(csd,file="data/fake_replicate2_csdata.RData")
  # csd=generate_fake_dataset(signal=F,eC=-4,biases.ref=biases.ref,replicate="3",condition="WT")
  # save(csd,file="data/fake_replicate3_csdata.RData")
  # csd=generate_fake_dataset(signal=F,eC=-3.8,biases.ref=biases.ref,replicate="4",condition="WT")
  # save(csd,file="data/fake_replicate4_csdata.RData")
  # csd=generate_fake_dataset(signal=F,eC=-4.1,biases.ref=biases.ref,replicate="5",condition="WT")
  # save(csd,file="data/fake_replicate5_csdata.RData")
  csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.3,replicate="1",condition="KO")
  save(csd,file="data/fake_signal_replicate1_csdata.RData")
  # csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-6.1,replicate="2",condition="KO")
  # save(csd,file="data/fake_signal_replicate2_csdata.RData")
  # csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-3.7,replicate="3",condition="KO")
  # save(csd,file="data/fake_signal_replicate3_csdata.RData")
  # csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.2,replicate="4",condition="KO")
  # save(csd,file="data/fake_signal_replicate4_csdata.RData")
  # csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-5.2,replicate="5",condition="KO")
  # save(csd,file="data/fake_signal_replicate5_csdata.RData")

  #bin raw data
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

base.res=5000

cs <- normalize_binless(cs, restart=F, ncores=ncores, base.res=base.res, fit.signal = fit.signal)#, tol=1e-4)

save(cs,file=paste0("data/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_csnorm_optimized.RData"))

for (resolution in c(5000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  #cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  if (nbg+nsig>1) {
    cs=detect_binned_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
    cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores, ref=cs@experiments[1,name])
  }
  save(cs,file=paste0("data/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_csnorm_optimized.RData"))
}

if (F) {
  load(paste0("data/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_csnorm_optimized.RData"))
  
  #observed
  resolution=5000
  mat=get_binned_matrices(cs, resolution=resolution, group="all")
  plot_binless_matrix(mat,upper="observed",lower="observed")
  ggsave(filename=paste0("images/fake_",nbg,"bg_",nsig,"sig_observed.png"),width=7,height=5)
  
  #normalized
  mat=get_binned_matrices(cs, resolution=resolution, group="all")
  plot_binless_matrix(mat,upper="normalized+1",lower="normalized+1")
  ggsave(filename=paste0("images/fake_",nbg,"bg_",nsig,"sig_",fit.signal,"_binned_normalized.png"),width=7,height=5)
  
  load("data/fake_0bg_1sig_TRUE_csnorm_optimized.RData")
  csT=cs
  load("data/fake_0bg_1sig_FALSE_csnorm_optimized.RData")
  csF=cs
  matF=get_binned_matrices(csF, resolution=resolution, group="all")
  matT=get_binned_matrices(csT, resolution=resolution, group="all")
  ice=merge(matF[,.(name,begin1,begin2,observed)],iterative_normalization(matF,niterations = 10),by=c("name","begin1","begin2"))
  mat=rbind(matF[,.(model="raw",begin1,begin2,value=observed)],
            matF[,.(model="binless w/o signal",begin1,begin2,value=normalized)],
            matT[,.(model="binless with signal",begin1,begin2,value=normalized)],
            ice[,.(model="equal coverage",begin1,begin2,value=ice.10)])
  mat[,model:=ordered(model,c("raw","equal coverage","binless w/o signal","binless with signal"))]
  mat[,avg:=.SD[begin1==0&begin2>0,mean(value)],by=model]
  mat[,value:=value/avg]
  
  #normalized side by side
  ggplot(mat)+
    geom_raster(aes(begin1,begin2,fill=value))+
    geom_raster(aes(begin2,begin1,fill=value))+coord_fixed()+facet_grid(~model)+
    scale_fill_gradient(low="white",high="black",na.value="white",trans="log10")+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme_minimal()+ theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"))
  ggsave(filename=paste0("images/fake_0bg_1sig_sidebyside_binned_normalized.pdf"),width=10,height=5)
  
  #column sums side by side
  colsums=rbind(mat[,.(model,begin1,begin2,value)],mat[begin1!=begin2,.(model,begin1=begin2,begin2=begin1,value)])[
    ,.(coverage=sum(value)),by=c("model","begin1")]
  ggplot(colsums)+geom_line(aes(begin1,coverage))+facet_grid(~model)+theme_minimal()+
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
    theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.spacing=unit(0,"cm"), panel.grid=element_blank())
  ggsave(filename=paste0("images/fake_0bg_1sig_sidebyside_binned_coverage.pdf"),width=10,height=2)
  
  
}

