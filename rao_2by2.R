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
ncores=4
base.res=5000

setwd("/home/yannick/simulations/cs_norm")

fit=T
bin=T
plot=T

if (fit==T) {
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
  save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
}

if (bin==T) {
  foreach (resolution=c(5000,10000),.errorhandling = "pass") %do% {
    cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
    cs=detect_binned_interactions(cs, resolution=resolution, group="all", ncores=ncores)
    cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=ncores)
    ref=cs@experiments[1,name]
    cs=detect_binned_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
    cs=detect_binless_differences(cs, ref, resolution=resolution, group="all", ncores=ncores)
  }
  #foreach (resolution=c(5000,10000),.errorhandling = "pass") %do% {
  #  cs=group_datasets(cs, resolution=resolution, verbose=T, ncores=ncores, group="condition")
  #  cs=detect_binned_interactions(cs, resolution=resolution, group="condition", ncores=ncores)
  #  cs=detect_binless_interactions(cs, resolution=resolution, group="condition", ncores=ncores)
  #  ref=get_binned_matrices(cs,resolution,group="condition")[,name[1]]
  #  cs=detect_binned_differences(cs, ref, resolution=resolution, group="condition", ncores=ncores)
  #  cs=detect_binless_differences(cs, ref, resolution=resolution, group="condition", ncores=ncores)
  #}
  save(cs,file=paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
}

if (plot==T) {
  #load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k_minlambda25.RData"))
  resolution=5000
  mat=get_binless_interactions(cs, resolution)
  plot_binless_matrix(mat,upper="log(observed)",lower="log(observed)")
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_observed.pdf"),width=15,height=15)
  plot_binless_matrix(mat,upper="log(binless)",lower="log(binless)")
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless.pdf"),width=15,height=15)
  plot_binless_signal_matrix(mat)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_signal.pdf"),width=15,height=15)
  mat=get_binless_differences(cs, cs@experiments[1,name],resolution)
  plot_binless_difference_matrix(mat)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_difference.pdf"),width=15,height=8)
  #
  #mat=get_binless_interactions(cs, resolution, group="condition")
  #plot_binless_matrix(mat,upper="log(observed)",lower="log(observed)")
  #ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_observed_grp.pdf"),width=10,height=10)
  #plot_binless_matrix(mat,upper="log(binless)",lower="log(binless)")
  #ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless_grp.pdf"),width=10,height=10)
  #plot_binless_signal_matrix(mat)
  #ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_signal_grp.pdf"),width=10,height=10)
  #ref=get_binned_matrices(cs,resolution,group="condition")[,name[1]]
  #mat=get_binless_differences(cs, ref, resolution, group="condition")
  #plot_binless_difference_matrix(mat)
  #ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_difference_grp.pdf"),width=10,height=10)
  
}

if (F) {
  #load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base",base.res/1000,"k.RData"))
  #generate plots
  resolution=5000
  
  #side-by-side at 5k: observed
  mat=get_binned_matrices(cs, resolution=resolution, group="all")
  plot_binless_matrix(mat,upper="observed",lower="observed")
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_observed.png"),width=10,height=10)
  
  #side-by-side at 5k: signal
  mat=get_binned_interactions(cs, resolution=resolution, group="all")
  mat[,is.significant:=prob.gt.expected>1-1e-12]
  plot_binless_signal_matrix(mat) + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_signal.png"),width=10,height=10)
  
  #side-by-side at 5k: normalized
  mat=get_binned_interactions(cs, resolution=resolution, group="all")
  mat[,is.significant:=prob.gt.expected>1-1e-12]
  plot_binless_matrix(mat, upper="normalized", lower="normalized") + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_normalized.png"),width=10,height=10)
  
  #side-by-side at 5k: differences
  mat=get_binned_differences(cs, resolution=resolution, group="all", ref=cs@experiments[1,name])
  mat[,is.significant:=pmax(`prob.gt.TBX3 GM12878 odd`,1-`prob.gt.TBX3 GM12878 odd`)>1-1e-12]
  mat[,direction:=ordered(direction,c("enriched","depleted"))]
  plot_binless_difference_matrix(mat) + geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_differences.png"),width=15,height=5)
  
  #side-by-side at 5k: binless signal
  mat=get_binless_interactions(cs)
  plot_binless_signal_matrix(mat)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless_signal.png"),width=10,height=10)
  
  #side-by-side at 5k: binless matrix
  mat=get_binless_interactions(cs)
  plot_binless_matrix(mat, upper="binless", lower="binless")
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binless.png"),width=10,height=10)
  
  #side-by-side at 5k: binless difference
  mat=get_binless_differences(cs,ref=cs@experiments[1,name])
  plot_binless_difference_matrix(mat)
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base10at",resolution/1000,"_binless_differences.png"),width=15,height=5)

  
  resolution=10000
  
  #side-by-side at 10k: observed
  mat=get_binned_matrices(cs, resolution=resolution, group="all")
  plot_binless_matrix(mat,upper="observed",lower="observed")
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_observed.png"),width=10,height=10)
  
  #side-by-side at 10k: signal
  mat=get_binned_interactions(cs, resolution=resolution, group="all")
  mat[,is.significant:=prob.gt.expected>1-1e-12]
  plot_binless_signal_matrix(mat) + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_signal.png"),width=10,height=10)
  
  #side-by-side at 10k: normalized
  mat=get_binned_interactions(cs, resolution=resolution, group="all")
  mat[,is.significant:=prob.gt.expected>1-1e-12]
  plot_binless_matrix(mat, upper="normalized", lower="normalized") + geom_point(aes(begin2,begin1),colour=muted("yellow"),data=mat[is.significant==T])
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_normalized.png"),width=10,height=10)
  
  #side-by-side at 10k: differences
  mat=get_binned_differences(cs, resolution=resolution, group="all", ref=cs@experiments[1,name])
  mat[,is.significant:=pmax(`prob.gt.TBX3 GM12878 odd`,1-`prob.gt.TBX3 GM12878 odd`)>1-1e-12]
  mat[,direction:=ordered(direction,c("enriched","depleted"))]
  plot_binless_difference_matrix(mat) + geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])
  ggsave(filename=paste0("images/rao_2by2_",sub,"_base5at",resolution/1000,"_binned_differences.png"),width=15,height=5)
  
  
  #get statistics on signal matrices
  resolution=5000
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
    data=rbind(data,
               dt[,.(name=sub,type="binned",resolution=resolution,variable=c("NGM","N1","N2","NIMR"),value=inter)],
               dt[,.(name=sub,type="binned",resolution=resolution,variable=c("JGM","J1","J2","JIMR"),value=jacard)])
    
    #jacard for binless interactions
    idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
    csg=cs@groups[[idx1]]
    idx2=get_cs_interaction_idx(csg, type="CSbsig", raise=T)
    csi=csg@interactions[[idx2]]
    csi@settings$min.l10FC=0.2
    mat=get_interactions(cs, type="CSbsig", resolution=resolution, group="all")
    mat[,value:=phi]
    mat = binless:::detect_binless_patches(mat, csi@settings)
    mat[,value:=NULL]
    #number of significant interactions, and percentage
    dt=mat[,.(num=sum(is.maximum)),by=name][
      ,.(name=locus,type="binless",resolution=resolution,variable=c("GM1","GM2","IMR1","IMR2"),value=num)]
    data=rbind(data,dt)
    #number of significant interactions
    compnames=data.table(t(matrix(cs@experiments[,name][c(c(1,2),c(1,3),c(2,4),c(3,4))],nrow=2)))
    dt = foreach (nref=compnames[,V1],n=compnames[,V2],.combine=rbind) %do% {
      refmat=mat[name==nref,.(bin1,bin2,ref.is=is.maximum)]
      merge(mat[name==n,.(name,bin1,bin2,is=is.maximum)],refmat,by=c("bin1","bin2"))[
        ,.(ref=nref,N=sum(is),refN=sum(ref.is),inter=sum(is==T&ref.is==T),union=sum(is==T|ref.is==T),
           jacard=sum(is==T&ref.is==T)/sum(is==T|ref.is==T)),by=name]
    }
    data=rbind(data,
               dt[,.(name=locus,type="binless",resolution=resolution,variable=c("NGM","N1","N2","NIMR"),value=inter)],
               dt[,.(name=locus,type="binless",resolution=resolution,variable=c("JGM","J1","J2","JIMR"),value=jacard)])
    dcast(data,name+type+resolution~variable,value.var = "value")


    
    
      dcast(melt(info[type=="binned"],id.vars = c("name","type","resolution"))[
    ,mean(value,na.rm=T),by=c("type","resolution","variable")],type+resolution~variable,value.var="V1")
  #enrichment vs detection type
  enrich=info[,.(name,resolution,type,Jdset=c(JGM,JIMR),Jrep=c(J1,J2))][(Jdset+Jrep)>0]
  enrich[,.N,by=c("resolution","type")]
  ggplot(enrich)+geom_boxplot(aes(type,Jdset/Jrep,colour=factor(resolution)))+geom_jitter(aes(type,Jdset/Jrep,colour=factor(resolution)))
  #jacard values at 5kb
  ggplot(melt(info[,.(name,type,resolution,JGM+JIMR,J1+J2)],id.vars=c("name","type","resolution"))[!is.na(value)&value>0])+
           geom_boxplot(aes(type,value,colour=factor(resolution)))+facet_wrap(~variable)
  
  
  
  
  #get statistics on difference matrices
  info = foreach(sub=subs, .combine=rbind, .errorhandling="remove") %dopar% {
    locus=strsplit(sub,"_")[[1]][1]
    load(paste0("data/rao_HiC_2by2_",sub,"_csnorm_optimized_base10k_dfuse",5,"_qmin_",0.01,"_stripped.RData"))
    foreach (resolution=c(5000,10000),.combine=rbind , .errorhandling="remove") %do% {
      data=data.table()
      
      #count number of interactions: binned GM vs GM
      mat=get_interactions(cs, resolution=resolution, group="all", type="CSdiff", ref=cs@experiments[1,name])
      mat=mat[name==cs@experiments[2,name]][(!is.na(is.significant))]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(is.significant),n.total=.N),by=name][
        ,.(name=locus,type="binned",resolution=resolution,cell="GM",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
      #count number of interactions: binned IMR vs IMR
      mat=get_interactions(cs, resolution=resolution, group="all", type="CSdiff", ref=cs@experiments[3,name])
      mat=mat[name==cs@experiments[4,name]][(!is.na(is.significant))]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(is.significant),n.total=.N),by=name][
        ,.(name=locus,type="binned",resolution=resolution,cell="IMR",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
      #count number of interactions: binnless GM vs GM
      idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
      csg=cs@groups[[idx1]]
      idx2=get_cs_interaction_idx(csg, type="CSbdiff", raise=T, ref=cs@experiments[1,name])
      csi=csg@interactions[[idx2]]
      csi@settings$min.l10FC=0.2
      mat=get_interactions(cs, type="CSbdiff", resolution=resolution, group="all", ref=cs@experiments[1,name])
      mat[,value:=delta]
      mat = binless:::detect_binless_patches(mat, csi@settings)
      mat[,value:=NULL]
      mat=mat[name==cs@experiments[2,name]]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(is.maximum||is.minimum),n.total=.N),by=name][
        ,.(name=locus,type="binless",resolution=resolution,cell="GM",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
      #count number of interactions: binnless GM vs GM
      idx1=get_cs_group_idx(cs, resolution, "all", raise=T)
      csg=cs@groups[[idx1]]
      idx2=get_cs_interaction_idx(csg, type="CSbdiff", raise=T, ref=cs@experiments[3,name])
      csi=csg@interactions[[idx2]]
      csi@settings$min.l10FC=0.2
      mat=get_interactions(cs, type="CSbdiff", resolution=resolution, group="all", ref=cs@experiments[3,name])
      mat[,value:=delta]
      mat = binless:::detect_binless_patches(mat, csi@settings)
      mat[,value:=NULL]
      mat=mat[name==cs@experiments[4,name]]
      #number of significant interactions
      dt=mat[,.(n.signif=sum(is.maximum||is.minimum),n.total=.N),by=name][
        ,.(name=locus,type="binless",resolution=resolution,cell="IMR",variable=c("N","pc"),value=c(n.signif,100*n.signif/n.total))]
      data=rbind(data,dt)
      
    }
  }
  info=info[name!="Comparison"]
  ggplot(info[variable=="pc"])+geom_boxplot(aes(type,value,colour=cell))+facet_wrap(~resolution)+scale_y_log10()
  ggsave(filename=paste0("images/rao_plot_FISH_signal_comparison_5k.png"),width=6,height=5)
  
  info[variable=="pc",.(num=.N,min=min(value),mean=mean(value),max=max(value)),by=c("resolution","type")]
  info[variable=="N",.(num=.N,min=min(value),mean=mean(value),max=max(value)),by=c("resolution","type")]
  
}

