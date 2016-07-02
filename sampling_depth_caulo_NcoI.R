library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")

### read caulobacter dataset and generate with different sampling depths

a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                  skip="SRR",nrows=1000000)
csd=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                     "data/caulo_NcoI_all10M", "WT", "1", skip="SRR", circularize=4042929, dangling.L=c(0,3,5), nrows=10000000,
                     dangling.R=c(3,0,-2), maxlen=600, save.data=F)
cs=merge_cs_norm_datasets(list(csd))
save(cs, file="data/caulo_NcoI_all10M_csnorm.RData")



### normalize different datasets

load("data/caulo_NcoI_all_csnorm.RData")

coverage=4
square.size=150000
bf_per_kb=0.25
cs=run_split_parallel(cs, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                      bf_per_decade=5, distance_bins_per_decade=100, verbose = F, iter=10000, ncores=30,
                      homogenize=F, outprefix="tmp/test")#, ops.count=ops.count, ops.bias=ops.bias)
#cs=run_split_parallel_recovery(cs, "tmp/test", square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
#                      bf_per_decade=5, distance_bins_per_decade=100, verbose = F, iter=10000, ncores=30,
#                      homogenize=F)
cs=postprocess(cs, resolution=10000, ncores=30, verbose=F)
cs@binned[[1]]=iterative_normalization(cs@binned[[1]], niterations=1)
save(cs, file="data/caulo_NcoI_all_csnorm_optimized.RData")
#save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))



### generate plots
prefix="NcoI_all"
samplings=c("1M","5M","10M")

#lFC
lFC = foreach (i=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=i,lFC)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
lFC=rbind(lFC,get_cs_binned(cs,1,"CS")[,.(dset="all",lFC)])
lFC[,dset:=ordered(dset,levels=c("all",samplings))]
ggplot(lFC)+geom_density(aes(lFC,colour=dset))
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_lFC.png"), width=10, height=7.5)

#normalized matrices
mat = foreach (i=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  get_cs_binned(cs,1,"CS")[,.(dset=i,begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)]
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
mat=rbind(mat,get_cs_binned(cs,1,"CS")[,.(dset="all",begin1,begin2,normalized,is.interaction,prob.observed.gt.expected)])
mat[,dset:=ordered(dset,levels=c("all",samplings))]
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")+
  facet_wrap(~dset)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_normalized.png"), width=10, height=7.5)

#nu
nu = foreach (i=samplings,.combine=rbind) %do% {
  load(paste0("data/caulo_",prefix,i,"_csnorm_optimized.RData"))
  data.table(pos=cs@biases[,pos],log_nu=cs@par$log_nu,dset=i)
}
load(paste0("data/caulo_",prefix,"_csnorm_optimized.RData"))
nu=rbind(nu,data.table(pos=cs@biases[,pos],log_nu=cs@par$log_nu,dset="all"))
nu[,dset:=ordered(dset,levels=c("all",samplings))]
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=dset))+xlim(1e6,2e6)
ggsave(filename=paste0("images/caulo_",prefix,"_sampling_depth_nu.png"), width=10, height=7.5)


### double run plots
#genomic biases
ggplot(melt(data.table(id=cs@biases[,pos],old=cs2@par$log_nu,new=cs@par$log_nu),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu.png"), width=10, height=7.5)
ggplot(melt(data.table(id=cs@biases[,pos],old=cs2@par$log_delta,new=cs@par$log_delta),id.vars="id",variable.name="delta"))+geom_line(aes(id,value,colour=delta))
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_delta.png"), width=10, height=7.5)
#ggplot(rbind(data.table(bias="nu",serial=opserial$par$log_nu,parallel=oppar$par$log_nu),data.table(bias="delta",serial=opserial$par$log_delta,parallel=oppar$par$log_delta)))+
#  geom_point(aes(serial,parallel))+facet_grid(.~bias)+stat_function(fun=identity)
#ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu_delta.png"), width=10, height=7.5)
#diagonal biases
ggplot(rbind(cs@par$decay[,.(dist,decay,dset="new")],cs2@par$decay[,.(dist,decay,dset="old")]))+
  geom_line(aes(dist,decay,colour=dset))+scale_x_log10()+scale_y_log10()
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_decay.png"), width=10, height=7.5)
c(cs2@par$eRJ,cs@par$eRJ)
c(cs2@par$eDE,cs@par$eDE)
c(cs2@par$eC,cs@par$eC)
c(cs2@binned[[1]]@alpha,cs@binned[[1]]@alpha)
cor(cs2@par$log_nu, cs@par$log_nu, method = "spearman")
cor(cs2@par$log_delta, cs@par$log_delta, method = "spearman")
cor(cs2@par$decay[,decay],cs@par$decay[,decay], method = "spearman")



### reload intermediate files of long run
#plot bias runs
ops=sapply(X=Sys.glob("ralph_*biases_op_*.RData"), FUN=function(x){ a=load(x); return(get(a[1]))}, USE.NAMES=T, simplify=F)
failures=data.table(failed=sapply(ops, function(x){length(grep("failed", tail(x$output,1)))>0}),
                    runtime=sapply(ops, function(x){x$runtime}),
                    sqsizes=tstrsplit(names(ops),"_")[[7]],
                    dispersion=sapply(ops, function(x){x$par$alpha}),
                    deviance=sapply(ops, function(x){x$par$deviance_proportion_explained}))
failures[,sqsizes:=as.numeric(sapply(sqsizes, function(x){substr(x,1,nchar(x)-6)}))]
ggplot(failures)+geom_point(aes(sqsizes, runtime, colour=failed))
ggplot(failures)+geom_point(aes(sqsizes, deviance, colour=failed))
ggplot(failures)+geom_point(aes(sqsizes, dispersion, colour=failed))+scale_y_log10()


#plot counts runs
ops=sapply(X=Sys.glob("tmp/ralph_*counts_op_*.RData"), FUN=function(x){ a=load(x); return(get(a[1]))}, USE.NAMES=T, simplify=F)
failures=data.table(failed=sapply(ops, function(x){length(grep("failed", tail(x$output,1)))>0}),
                    runtime=sapply(ops, function(x){x$runtime}),
                    sqsizes=tstrsplit(names(ops),"_")[[7]],
                    deviance=sapply(ops, function(x){x$par$deviance_proportion_explained}))
failures[,sqsizes:=as.numeric(sapply(sqsizes, function(x){substr(x,1,nchar(x)-6)}))]
failures[,idx:=.I]
ggplot(failures)+geom_point(aes(idx, runtime, colour=failed))
ggplot(failures)+geom_point(aes(idx, deviance, colour=failed))

ops.bias=sapply(X=Sys.glob(paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_biases_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
ops.bias = output.binder(ops.bias)
save(ops.bias, file=paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_ops_bias.RData"))

ops.count=sapply(X=Sys.glob(paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_counts_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
ops.count = output.binder(ops.count)
save(ops.count, file=paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_ops_count.RData"))


### generate LGF matrix
a=fread("~/simulations/caulobacter/data/GSM1120448_Laublab_NcoI_HiC_NA1000_swarmer_cell_untreated_overlap_after_HiCNorm_fullgen.count")
setnames(a,c("bin1","bin2","lgfcount"))
bins=data.table(begin=seq(1,4045000+10000,10000))
bins[,name:=paste0("bin",.I-1)]
setkey(bins,name)
setkey(a,bin1)
a=bins[a]
setnames(a,c("begin", "name"), c("begin1","bin1"))
setkey(a,bin2)
a=bins[a]
setnames(a,c("begin", "name"), c("begin2","bin2"))
#newbins=seq(1,4045000+20000,20000)
#a[,nb1:=cut2(begin1,newbins)]
#a[,nb2:=cut2(begin2,newbins)]
#a=a[,.(min(begin1),min(begin2),lgfcount=sum(lgfcount)),by=c("nb1","nb2")]
#a[,c("ign","begin1","ign2"):=tstrsplit(nb1,'[[,]')]
#a[,c("ign","begin2","ign2"):=tstrsplit(nb2,'[[,]')]
#a[,c("begin1","begin2"):=list(as.numeric(begin1),as.numeric(begin2))]
#ggplot(a[begin1<=2003200&begin2<=2003200])+geom_raster(aes(begin1,begin2,fill=log(lgfcount)))+geom_raster(aes(begin2,begin1,fill=log(lgfcount)))+
#  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")
ggplot(a)+geom_raster(aes(begin1,begin2,fill=log(lgfcount)))+geom_raster(aes(begin2,begin1,fill=log(lgfcount)))+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")
ggsave(filename = paste0("images/",prefix,"_HiCNorm.png"), width=10, height=9)+  theme(legend.position = "none") 


