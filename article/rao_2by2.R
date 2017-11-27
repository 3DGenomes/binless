library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

args=commandArgs(trailingOnly=TRUE)
sub=args[1]
ncores=1
base.res=5000

setwd("/home/yannick/simulations/cs_norm")

fit=T
restart=F
bin=T
plot=T

if (fit==T) {
  if (restart==F) {
    load(paste0("data/rao_HiCodd_GM12878_",sub,"_csdata.RData"))
    csd1=csd
    load(paste0("data/rao_HiCeven_GM12878_",sub,"_csdata.RData"))
    csd2=csd
    load(paste0("data/rao_HiCodd_IMR90_",sub,"_csdata.RData"))
    csd3=csd
    load(paste0("data/rao_HiCeven_IMR90_",sub,"_csdata.RData"))
    csd4=csd
    cs=merge_cs_norm_datasets(list(csd1,csd2,csd3,csd4), different.decays="none")
    cs <- normalize_binless(cs, restart=F, ncores=ncores, base.res=base.res)
  } else {
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
    if (!has_converged(cs)) cs <- normalize_binless(cs, restart=T, ngibbs = 15, ncores=ncores)
  }
  save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
}

if (bin==T) {
  load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  
  foreach (resolution=c(5000,10000),.errorhandling = "pass") %do% {
    cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
    cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
    cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
    ref=cs@experiments[1,name]
    cs=detect_binned_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
    cs=detect_binless_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
    ref=cs@experiments[3,name]
    cs=detect_binned_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
    cs=detect_binless_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
  }
  save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
}

if (plot==T) {
  load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  resolution=5000
  mat=get_binless_interactions(cs, resolution)
  #observed
  plot_binless_matrix(mat,upper="observed",lower="observed")+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_observed.pdf"),width=20,height=8)
  #binless and observed
  plot_binless_matrix(mat,upper="binless",lower="observed")+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless_observed.pdf"),width=20,height=8)
  #binless
  plot_binless_matrix(mat,upper="binless",lower="binless")+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless.pdf"),width=20,height=8)
  #binless signal
  plot_binless_signal_matrix(mat)+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless_signal.pdf"),width=20,height=8)
  #binless differences
  mat=get_binless_differences(cs, cs@experiments[1,name],resolution)
  plot_binless_difference_matrix(mat)+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless_differences.pdf"),width=16,height=8)
  #binned normalized
  mat=get_binned_interactions(cs, resolution=resolution, group="all")
  mat[,is.significant:=prob.gt.expected>=1-1e-12]
  plot_binless_matrix(mat, upper="normalized", lower="normalized") + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_normalized.pdf"),width=16,height=8)
  #binned differences
  mat=get_binned_differences(cs, resolution=resolution, group="all", ref=cs@experiments[1,name])
  mat[,direction:=ordered(direction,c("enriched","depleted"))]
  plot_binless_difference_matrix(mat) + geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+facet_wrap(~name,nrow = 1)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_differences.pdf"),width=16,height=8)
}

if (F) {
  
  subs=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "FOXP1_1.3M",
         "TBX3_1.5M", "Comparison_1.7M", "Talk_2M", "ADAMTS1_2.3M")
  registerDoParallel(cores=10)
  info = foreach(sub=subs,.combine=rbind,.errorhandling="remove") %dopar% {
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
    data.table(ori=sub,has.converged=binless:::has_converged(cs),#run=run,
               runtime=cs@diagnostics$params[,sum(runtime)],nsteps=cs@diagnostics$params[,max(step)],
               name=cs@experiments[,name],alpha=cs@par$alpha,
               lambda1=cs@par$lambda1,lambda2=cs@par$lambda2,eCprime=cs@par$eCprime,lambda_diag=cs@par$lambda_diag,
               ngroups=length(cs@groups))#,
  }
  info
  
  
  #get statistics on signal matrices
  resolution=10000
  data=data.table()
    
  #jacard for binned interactions
  mat=get_binned_interactions(cs, resolution=resolution, group="all")[!is.na(is.significant)]
  #number of significant interactions
  dt=mat[,.(num=sum(is.significant)),by=name][
    ,.(name=sub,type="binned",resolution=resolution,variable=c("GM1","GM2","IMR1","IMR2"),value=num)]
  data=rbind(data,dt)
  #number of significant interactions
  compnames=data.table(t(matrix(cs@experiments[,name][c(c(1,2),c(1,3),c(2,4),c(3,4))],nrow=2)))
  dt = foreach (nref=compnames[,V1],n=compnames[,V2],.combine=rbind) %do% {
    refmat=mat[name==nref,.(bin1,bin2,ref.is=is.significant)]
    merge(mat[name==n,.(name,bin1,bin2,is=is.significant)],refmat,by=c("bin1","bin2"))[
      ,.(ref=nref,N=sum(is),refN=sum(ref.is),inter=sum(is==T&ref.is==T),union=sum(is==T|ref.is==T),
         jacard=sum(is==T&ref.is==T)/sum(is==T|ref.is==T)),by=name]
  }
  
  data
  dt[,.(name=sub,type="binned",resolution=resolution,variable=c("NGM","N1","N2","NIMR"),value=inter)]
  dt[,.(name=sub,type="binned",resolution=resolution,variable=c("JGM","J1","J2","JIMR"),value=jacard)]
  
  #get statistics on difference matrices
  subs=c("SELP_150k","Peak1_450k","ADAMTS2_450k","PARM1_600k","Tbx19_700k","SEMA3C_1M", "FOXP1_1.3M",
         "TBX3_1.5M", "Comparison_1.7M", "Talk_2M", "ADAMTS1_2.3M")
  registerDoParallel(cores=10)
  info = foreach(sub=subs, .combine=rbind, .errorhandling="remove") %dopar% {
    locus=strsplit(sub,"_")[[1]][1]
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
    foreach (resolution=c(5000,10000),.combine=rbind , .errorhandling="remove") %do% {
      data=data.table()
      
      #count number of interactions: binned GM vs GM
      mat=get_binned_differences(cs, resolution=resolution, group="all", ref=cs@experiments[1,name])
      mat=mat[name==cs@experiments[2,name]][(!is.na(is.significant))]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(is.significant),n.total=.N),by=name][
        ,.(name=locus,type="binned",resolution=resolution,cell="GM",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
      #count number of interactions: binned IMR vs IMR
      mat=get_binned_differences(cs, resolution=resolution, group="all", ref=cs@experiments[3,name])
      mat=mat[name==cs@experiments[4,name]][(!is.na(is.significant))]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(is.significant),n.total=.N),by=name][
        ,.(name=locus,type="binned",resolution=resolution,cell="IMR",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
      #count number of interactions: binnless GM vs GM
      mat=get_binless_differences(cs, resolution=resolution, group="all", ref=cs@experiments[1,name])
      mat=mat[name==cs@experiments[2,name]]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(delta!=0),n.total=.N),by=name][
        ,.(name=locus,type="binless",resolution=resolution,cell="GM",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
      #count number of interactions: binnless GM vs GM
      mat=get_binless_differences(cs, resolution=resolution, group="all", ref=cs@experiments[3,name])
      mat=mat[name==cs@experiments[4,name]]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(delta!=0),n.total=.N),by=name][
        ,.(name=locus,type="binless",resolution=resolution,cell="IMR",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
    }
  }
  info=info[!(type=="binless"&resolution==10000)] #normalization was 5k base res, 10k inconsistent
  ggplot(info[type=="binned"&variable=="pc"])+geom_boxplot(aes(name,value))+geom_jitter(aes(name,value,colour=factor(resolution)))+scale_y_log10()
  ggplot(info[variable=="pc"])+geom_boxplot(aes(type,value,colour=cell))+facet_wrap(~resolution)+scale_y_log10()
  
  info[variable=="pc",.(num=.N,min=min(value),mean=mean(value),max=max(value)),by=c("resolution","type")]
  info[variable=="N",.(num=.N,min=min(value),mean=mean(value),max=max(value)),by=c("resolution","type")]
  
}

