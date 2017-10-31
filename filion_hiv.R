library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)


setwd("/home/yannick/simulations/cs_norm")

if(F) {
  a=examine_dataset("/scratch/filion/E11N_chr7_137-139M_imperfect.tsv",skip=0,nrows=1e6,skip.fbm=F, read.len=75)
  #cut at 750bp for pdiag
  #cut at 1000bp for pclose
  #dL=1 and dR=5
  csd=read_and_prepare(paste0("/scratch/filion/E11N_chr7_137-139M_imperfect.tsv"),
                       paste0("data/filion_E11N_chr7_1.6M"), "E11N", "1", locus=c("chr7",137200000,138800000),
                       enzyme="MboI", name=paste("E11N"), circularize=-1, dangling.L=c(1),
                       dangling.R=c(5), maxlen=750, read.len=75, dmin=1000, save.data=T)
  csd=read_and_prepare(paste0("/scratch/filion/E11P_chr7_137-139M_imperfect.tsv"),
                       paste0("data/filion_E11P_chr7_1.6M"), "E11P", "1", locus=c("chr7",137200000,138800000),
                       enzyme="MboI", name=paste("E11P"), circularize=-1, dangling.L=c(1),
                       dangling.R=c(5), maxlen=750, read.len=75, dmin=1000, save.data=T)
  
  load("data/filion_E11P_chr7_2M_csdata_with_data.RData")
  plot_raw(csd@data,b1=csd@data[,min(rbegin1)+2000],e1=csd@data[,min(rbegin1)+5000])
}

bpk=50
dfuse=20 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
ncores=10

setwd("/home/yannick/simulations/cs_norm")

load("data/filion_E11N_chr7_1.6M_csdata.RData")
csd1=csd
load("data/filion_E11P_chr7_1.6M_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none", dfuse=dfuse)
cs = normalize_binless(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 5, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
               ncounts = 1000000, ncores=ncores, base.res=10000, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
save(cs,file=paste0("data/filion_E11_chr7_1.6M_csnorm_optimized.RData"))
cs = normalize_binless(cs, restart=T, ngibbs = 15, ncores=ncores)
save(cs,file=paste0("data/filion_E11_chr7_1.6M_csnorm_optimized.RData"))

for (resolution in c(10000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  ref=cs@experiments[1,name]
  cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores)
  save(cs,file=paste0("data/filion_E11_chr7_1.6M_csnorm_optimized.RData"))
}

if (F) {
  load(paste0("data/filion_E11_chr7_1.6M_csnorm_optimized.RData"))
  
  binless:::has_converged(cs)
  cs@diagnostics$params[,sum(runtime)]/3600
  
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
  
  resolution=10000
  mat.binned=get_matrices(cs, resolution=resolution, group="all")
  #observed
  ggplot(mat.binned)+geom_raster(aes(begin1,begin2,fill=log(observed)))+geom_raster(aes(begin2,begin1,fill=log(observed)))+facet_wrap(~name)+
    scale_fill_gradient(high="black", low="white", na.value = "white")+coord_fixed()
  ggsave(filename="images/filion_E11_base10k_observed.pdf",width=20,height=7)
  
  mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi))+
    geom_raster(aes(begin2,begin1,fill=phi))+coord_fixed()+
    facet_wrap(~name)+scale_fill_gradient(high="black", low="white",na.value = "white")
  ggsave(filename="images/filion_E11_base10k_binless.pdf",width=20,height=7)
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=beta))+
    geom_raster(aes(begin2,begin1,fill=beta))+coord_fixed()+
    facet_wrap(~name)+scale_fill_gradient(high="black", low="white",na.value = "white")
  ggsave(filename="images/filion_E11_base10k_binless_nosignif.pdf",width=20,height=7)
  
  write.table(merge(mat.binned[,.(name,begin1,begin2,observed)],
                    mat[,.(name,begin1,begin2,phi,beta)],by=c("name","begin1","begin2")),
              file="data/filion_E11_base10k_normalized.txt", row.names=F, quote=F)
  
  mat=get_interactions(cs, type="CSbdiff", resolution=resolution, group="all", ref=as.character(cs@experiments[1,name]))
  ggplot(mat)+geom_raster(aes(begin1,begin2,fill=delta))+coord_fixed()+
    geom_raster(aes(begin2,begin1,fill=delta))+
    facet_wrap(~name)+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"),na.value="white")
  ggsave(filename="images/filion_E11_base10k_bdiff.pdf",width=15,height=7)
  
}
