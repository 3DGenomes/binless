library(binless)
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
  mat = foreach(name=c("FISH1","FISH2","FISH3","FISH4"), size=c("1.5M","1.5M","1.5M","2.2M")) %dopar% {
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
                  idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
                  csg=cs@groups[[idx1]]
                  idx2=get_cs_interaction_idx(csg, type="CSbsig", raise=T)
                  csi=csg@interactions[[idx2]]
                  csi@settings$min.l10FC=0.2
                  mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
                  mat[,value:=phi]
                  mat = binless:::detect_binless_patches(mat, csi@settings)
                  mat[,value:=NULL]
                  mat[,phi.max:=ifelse(is.maximum==T,NA,phi)]
                  mat=merge(mat,get_matrices(cs, resolution=resolution, group="all")[,.(name,bin1,bin2,decaymat)],by=c("name","bin1","bin2"))
                  p=ggplot(mat)+geom_raster(aes(begin1,begin2,fill=phi/log(10)+log10(decaymat)))+
                    geom_raster(aes(begin2,begin1,fill=phi.max/log(10)+log10(decaymat)))+
                    scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value="black")+
                    geom_segment(aes(x=V2,y=V2,xend=V3,yend=V3,colour=V4),data=locs,size=3)+
                    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) + guides(colour=F) +
                    theme_void()+ theme(axis.title=element_blank(),
                                        panel.background = element_rect(fill = "white", colour = "black"),
                                        panel.spacing=unit(0,"cm"))+
                    coord_fixed()+labs(fill="log10 FC")
                  ggsave(p,filename=paste0("images/rao_plot_",sub,"_",run,"_",resolution/1000,"k_binless_with_decay.png"),width=6,height=5)
            }
  }
  
  #get statistics
  run="cv_cvsd_outlier_rmdiag"
  info = foreach(name=c("FISH1","FISH2","FISH3","FISH4"), size=c("1.5M","1.5M","1.5M","2.2M"), .combine=rbind) %dopar% {
    sub=paste0(name,"_",size)
    load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base10k_bpk",bpk,"_dfuse",dfuse,"_cv_cvsd_outlier_rmdiag_stripped.RData"))
    locs=fread("/scratch/rao/mapped/interesting_locations_hg38_FISH_oligos.bed")[
      grep(paste0(name,"_L"),V4),.(Lbegin=V2,Lend=V3,loc=V4)]
    locs[,loc:=as.numeric(tstrsplit(loc,"_L")[[2]])]
    setkey(locs,Lbegin,Lend)
    foreach (resolution=c(5000,10000), .combine=rbind) %do% {
      mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
      mat=merge(mat,get_interactions(cs, resolution=resolution, group="all",type="CSsig")[
        ,.(name,bin1,bin2,signal,is.significant,prob.gt.expected)],by=c("name","bin1","bin2"))
      mat[,name:=NULL]
      mat=foverlaps(mat,locs,by.x=c("begin1","end1"), type="any")
      setnames(mat,c("Lbegin","Lend","loc"),c("Lbegin1","Lend1","loc1"))
      mat=foverlaps(mat,locs,by.x=c("begin2","end2"), type="any")
      setnames(mat,c("Lbegin","Lend","loc"),c("Lbegin2","Lend2","loc2"))
      dt1=mat[loc1+loc2==3,.(name=name,resolution=resolution,locus="L1+L2 (positive)",is.maximum,phi,patchno,signal,is.significant,prob.gt.expected)]
      dt2=mat[loc1+loc2==5,.(name=name,resolution=resolution,locus="L2+L3 (negative)",is.maximum,phi,patchno,signal,is.significant,prob.gt.expected)]
      rbind(dt1,dt2)
    }
  }
  
  #binless
  ggplot(info)+geom_boxplot(aes(locus,phi,colour=factor(resolution)))+facet_wrap(~name)
  #binned
  ggplot(info)+geom_boxplot(aes(locus,log(signal),colour=factor(resolution)))+facet_wrap(~name)
  
  info.vs=melt(info[resolution==5000&signal>0,.(locus,binned=log10(signal),binless=phi/log(10),name)],
               measure.vars = c("binned","binless"), variable.name="method", value.name="log10 FC")
  ggplot(info.vs)+geom_boxplot(aes(locus,`log10 FC`,colour=`method`))+facet_wrap(~name)
  ggsave(filename=paste0("images/rao_plot_FISH_signal_comparison_5k.png"),width=6,height=5)
  
  info.vs[,.N,by=c("locus","method","name")]
        
  
  #binless
  ggplot(info[,.(nsignif=sum(is.maximum)),by=c("locus","resolution","name")])+
    geom_point(aes(locus,nsignif,colour=factor(resolution)))+facet_wrap(~name)
  #binned
  ggplot(info[,.(nsignif=sum(is.significant)),by=c("locus","resolution","name")])+
    geom_point(aes(locus,nsignif,colour=factor(resolution)))+facet_wrap(~name)
  ggplot(info)+geom_boxplot(aes(locus,prob.gt.expected,colour=factor(resolution)))+facet_wrap(~name)
  
