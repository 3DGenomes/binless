library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)

args=commandArgs(trailingOnly=TRUE)
sub=args[1]
bpk=30
dfuse=5 #as.integer(args[2])
bpd=10 #as.integer(args[4])
bpb=10 #as.integer(args[5])
ncores=5
qmin=0.01

setwd("/home/yannick/simulations/cs_norm")

if (FALSE) {
  a=examine_dataset("/scratch/rao/mapped/GM12878_MboI_in_situ/GM12878_MboI_HIC003_FISH1.tsv", skip=0L,nrows=1000000, skip.fbm=F,
                    read.len=101)
  foreach(name=c("FISH1","FISH2","FISH3","FISH4"), size=c("1.5M","1.5M","1.5M","2.2M")) %do% {
   csd=read_and_prepare(paste0("/scratch/rao/mapped/GM12878_MboI_in_situ/GM12878_MboI_HICall_",name,".tsv"),
                        paste0("data/rao_HiCall_GM12878_",name,"_",size), "GM", "1",
                        enzyme="MboI", name=paste(name,"GM12878 all"), circularize=-1, dangling.L=c(0),
                        dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)
  }
  load("data/rao_HiCall_GM12878_FISH4_2.2M_csdata_with_data.RData")
  data=get_raw_reads(csd@data, csd@biases[,min(pos)], csd@biases[,max(pos)])
  plot_binned(data, resolution=10000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])
  #plot_raw(data[rbegin2<min(rbegin1)+10000])
  
}


#load(paste0("data/rao_HiCall_GM12878_",sub,"_csdata.RData"))
#csd1=csd
#cs=merge_cs_norm_datasets(list(csd1), different.decays="none", dfuse=dfuse, qmin=qmin)
#cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
#               ngibbs = 25, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-2, tol.leg=1e-4,
#               ncounts = 1000000, ncores=ncores, base.res=10000, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)
#save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag.RData"))

load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag.RData"))

for (resolution in c(10000)) {
  cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
  cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
  save(cs,file=paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag.RData"))
}

if(F) {
  #generate plots for all datasets
  run="cv_cvsd_outlier_rmdiag"
  mat = foreach(name=c("FISH1","FISH2","FISH3","FISH4"), size=c("1.5M","1.5M","1.5M","2.2M")) %do% {
            sub=paste0(name,"_",size)
            load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag_stripped.RData"))
            locs=fread("/scratch/rao/mapped/interesting_locations_hg38_FISH_oligos.bed")[grep(paste0(name,"_L"),V4)]
            foreach (resolution=c(5000,10000)) %do% {
                  #observed
                  mat.raw=get_matrices(cs, resolution=resolution, group="all")
                  p=ggplot(mat.raw)+
                    geom_raster(aes(begin1,begin2,fill=observed))+
                    geom_raster(aes(begin2,begin1,fill=observed))+coord_fixed()+
                    scale_fill_gradient(low="white",high="black",na.value="white",trans="log")+
                    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
                    theme_minimal()+ theme(axis.title=element_blank(),
                                           panel.background = element_rect(fill = "white", colour = "black"),
                                           panel.spacing=unit(0,"cm"))
                  ggsave(p,filename=paste0("images/rao_plot_",sub,"_",run,"_",resolution/1000,"k_observed.png"),width=6,height=5)
                  
                  #binless with decay
                  mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
                  mat=merge(mat,mat.raw[,.(name,bin1,bin2,decaymat)],by=c("name","bin1","bin2"))
                  a=mat[is.maximum==T]
                  a=a[,.SD[,.(begin1=c(begin1,begin1,end1,end1)-resolution/2, begin2=c(begin2,end2,begin2,end2)-resolution/2,
                              patchno, phi)][chull(begin1,begin2)], by=c("patchno","name")]
                  p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi/log(10)+log10(decaymat)))+
                    geom_raster(aes(begin2,begin1,fill=phi/log(10)+log10(decaymat)))+
                    geom_polygon(aes(begin2,begin1,group=patchno),colour="black",fill=NA,data=a)+
                    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="white")+
                    geom_segment(aes(x=V2,y=V2,xend=V3,yend=V3,colour=V4),data=locs,size=3)+
                    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
                    theme_minimal()+ theme(axis.title=element_blank(),
                                        panel.background = element_rect(fill = "white", colour = "black"),
                                        panel.spacing=unit(0,"cm"))+coord_fixed()+labs(fill="log10 FC")
                  ggsave(p,filename=paste0("images/rao_plot_",sub,"_",run,"_",resolution/1000,"k_binless_with_decay.png"),width=6,height=5)
            }
  }
}
