library(ggplot2)
library(data.table)
library(csnorm)

setwd("/home/yannick/simulations/cs_norm")


#read_and_prepare("/scratch/ledily/mapped/HindIII_T0_both_filled_map_chr6.tsv", "data/ledily_T0_HindIII_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))
#read_and_prepare("/scratch/ledily/mapped/HindIII_T60_both_filled_map_chr6.tsv", "data/ledily_T60_HindIII_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))

#read_and_prepare("/scratch/ledily/mapped/NcoI_T0_both_filled_map_chr6.tsv", "data/ledily_T0_NcoI_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))
#read_and_prepare("/scratch/ledily/mapped/NcoI_T60_both_filled_map_chr6.tsv", "data/ledily_T60_NcoI_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))

#read_and_prepare("/scratch/rao/mapped/HICall_both_filled_map_chr1.tsv", "data/rao_HICall_chr1_all", skip="SRR")

#read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
#                 "data/caulo_BglIIr1_all", skip="SRR", circularize=4042929)
#read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate2_reads_int.tsv",
#                 "data/caulo_BglIIr2_all", skip="SRR", circularize=4042929)
a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                  skip="SRR",nrows=1000000)
csd=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                 "data/caulo_NcoI_all", "WT", "1", skip="SRR", circularize=4042929, dangling.L=c(0,3,5),
                 dangling.R=c(3,0,-2), maxlen=600, save.data=F)
cs=merge_cs_norm_datasets(list(csd))

#read_and_prepare("/scratch/ralph/HiC/3_Mapped/Bcell_Sox2_10Mb_both_filled_map.tsv", "data/ralph_Bcell_Sox2", skip="HWI")
#read_and_prepare("/scratch/ralph/HiC/3_Mapped/EScell_Sox2_10Mb_both_filled_map.tsv", "data/ralph_EScell_Sox2", skip="HWI")

dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = 0, eRJ = 5, eDE = 7)
counts=dset$counts
biases=dset$biases
dset_statistics(biases,counts)

#prefix="ledily_T60_HindIII_chr6"
#prefix="caulo_NcoI_all"
prefix="caulo_NcoI_1004680-2998140"
#prefix="ralph_Bcell_Sox2"
biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))

#ledily:  tad 821 151000000-153000000
biases=biases[pos>=151000000&pos<=153000000]
#biases=biases[pos>=150000000&pos<=154000000]

#ralph: restrict to 30500000,32500000
#biases=biases[pos<=30600714]
#biases=biases[pos>=30500000&pos<=32500000]
#biases=biases[pos>=30500000&pos<=30600000]
beginrange=biases[1,id]
endrange=biases[.N,id]
biases[,id:=id-beginrange+1]
prefix=biases[,paste0(prefix,"_",min(pos),"-",max(pos))]
counts=counts[id1>=beginrange&id1<=endrange&id2>=beginrange&id2<=endrange]
counts[,c("id1","id2"):=list(id1-beginrange+1,id2-beginrange+1)]

#caulo/toy: restrict to 100 consecutive rsites
biases=biases[pos>=1000000&pos<=3000000]
beginrange=biases[1,id]
endrange=biases[.N,id]
beginrange=1
endrange=biases[pos<4042929/2][,.N]
beginrange=1
endrange=975
beginrange=1
endrange=200
beginrange=1
endrange=80
beginrange=22
endrange=96
beginrange=49
endrange=119
biases=biases[beginrange:endrange]
biases[,id:=id-beginrange+1]
prefix=biases[,paste0("caulo_NcoI_",min(pos),"-",max(pos))]
counts=counts[id1>=beginrange&id1<=endrange&id2>=beginrange&id2<=endrange]
counts[,c("id1","id2"):=list(id1-beginrange+1,id2-beginrange+1)]


counts=fill_zeros(counts,biases)
counts[,distance:=abs(pos2-pos1)] #not filled
write.table(biases, file = paste0("data/",prefix,"_biases.dat"), quote=F, row.names = F)
write.table(counts, file = paste0("data/",prefix,"_counts.dat"), quote=F, row.names = F)

counts[,distance:=pmin(abs(pos2-pos1),4042929-abs(pos2-pos1)+1)] #not filled


### effect of parallelization
opall <- csnorm_fit(biases, counts, bf_per_kb=1,
                                bf_per_decade=5, verbose = T, iter=100000)
opall=postprocess(biases, counts, opall, resolution=10000, ncores=30)
#save(opall, file = paste0("data/",prefix,"_op_maxcount_-1_redo.RData"))
load(paste0("data/",prefix,"_op_maxcount_-1.RData"), verbose=T)

coverage=4
square.size=150000
oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=1,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = F, iter=100000, ncores=30, homogenize=F)
oppar=postprocess(biases, counts, oppar, resolution=10000, ncores=30)
#save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb10.RData"))

load("data/",prefix,"_op_maxcount_-1.RData")
opserial=op
load(paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k.RData"))

#lFC histogram and matrices
ggplot(rbind(opserial$mat[,.(bin1,bin2,dset="serial",lFC)],oppar$mat[,.(bin1,bin2,dset="parallel",lFC)]))+geom_histogram(aes(lFC,fill=dset),position="dodge")
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_hist.png"), width=10, height=7.5)
ggplot()+geom_raster(data=oppar$mat,aes(begin1,begin2,fill=lFC))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+
  geom_raster(data=opserial$mat,aes(begin2,begin1,fill=lFC))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opserial$mat[is.interaction==T])
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_mat.png"), width=10, height=7.5)
#normalized matrix
ggplot()+geom_raster(data=oppar$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+
  geom_raster(data=opserial$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opserial$mat[is.interaction==T])
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_normalized.png"), width=10, height=7.5)
#genomic biases
#ggplot(melt(data.table(id=biases[,pos],serial=exp(opserial$par$eRJ+opserial$par$log_nu),parallel=exp(oppar$par$eRJ+oppar$par$log_nu)),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))+geom_point(data=biases,aes(pos,rejoined))
ggplot(melt(data.table(id=biases[,pos],serial=opserial$par$log_nu,parallel=oppar$par$log_nu),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu.png"), width=10, height=7.5)
ggplot(melt(data.table(id=biases[,pos],serial=opserial$par$log_delta,parallel=oppar$par$log_delta),id.vars="id",variable.name="delta"))+geom_line(aes(id,value,colour=delta))
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_delta.png"), width=10, height=7.5)
ggplot(rbind(data.table(bias="nu",serial=opserial$par$log_nu,parallel=oppar$par$log_nu),data.table(bias="delta",serial=opserial$par$log_delta,parallel=oppar$par$log_delta)))+
  geom_point(aes(serial,parallel))+facet_grid(.~bias)+stat_function(fun=identity)
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu_delta.png"), width=10, height=7.5)
#diagonal biases
ggplot(melt(data.table(dist=counts[,distance],serial=opserial$pred$log_decay_close,parallel=oppar$pred$log_decay_far, key="dist"),id.vars="dist"))+
  geom_line(aes(dist,value,colour=variable))+scale_x_log10()
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_decay.png"), width=10, height=7.5)
c(opserial$mat[,mean(lFC)], oppar$mat[,mean(lFC)])
c(opserial$par$eRJ,oppar$par$eRJ)
c(opserial$par$eDE,oppar$par$eDE)
c(opserial$par$eC,oppar$par$eC)
c(opserial$disp$alpha,oppar$disp$alpha)




### single run plots
coverage=4
square.size=150000
bf_per_kb=0.25
cs=run_split_parallel(cs, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = F, iter=10000, ncores=30,
                         homogenize=F, outprefix="tmp/test")
#oppar=run_split_parallel_recovery(counts, biases, outprefix, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
#                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=10000, ncores=30, homogenize=F, circularize=circularize)
cs=postprocess(cs, resolution=10000, ncores=30, verbose=F)
oppar$ice=iterative_normalization(oppar$mat, niterations=1)
save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))

load("data/caulo_NcoI_3189-2020271_op_maxcount_-1_parallel_inhomogeneous_cov4X_sq150k_bfpkb01.RData")
load("data/caulo_BglIIr1_1-2009521_op_maxcount_-1_parallel_inhomogeneous_cov4X_sq150k_bfpkb01.RData")
load("data/caulo_BglIIr2_15389-2009521_op_maxcount_-1_parallel_inhomogeneous_cov4X_sq150k_bfpkb01.RData")

#prefix="ledily_T0_NcoI_chr6_150000151-153992077"
#prefix="ledily_T0_NcoI_chr6_151001273-152989574"
#prefix="ledily_T0_HindIII_chr6_150004969-153999619"
#prefix="ledily_T0_HindIII_chr6_151000912-152997061"
#prefix="ledily_T60_NcoI_chr6_150000151-153992077"
#prefix="ledily_T60_NcoI_chr6_151001273-152989574"
#prefix="ledily_T60_HindIII_chr6_150004969-153999619"
#prefix="ledily_T60_HindIII_chr6_151000912-152997061"
#prefix="caulo_BglIIr1_1-2009521"
#prefix="caulo_BglIIr2_15389-2009521"
#prefix="caulo_NcoI_3189-2020271"
prefix="caulo_NcoI_all"
#prefix="caulo_NcoI_1004680-2998140"

biases=fread(paste0("data/",prefix,"_biases.dat"))
counts=fread(paste0("data/",prefix,"_counts.dat"))
load(paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))
setkey(oppar$mat, begin1,begin2)
#matrix and detected interactions
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(normalized)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none")
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_normalized.png"), width=10, height=9)
#lFC
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=lFC))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=lFC))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_lFC_mat.png"), width=10, height=9)
#raw
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(observed)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+scale_fill_gradient(low="white", high="black")+theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_observed.png"), width=10, height=9)
#ice
ggplot()+geom_raster(data=oppar$ice, aes(begin1,begin2,fill=log(N)))+geom_raster(data=oppar$ice, aes(begin2,begin1,fill=log(N)))+scale_fill_gradient(low="white", high="black")+  theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_ICE.png"), width=10, height=9)+  theme(legend.position = "none") 
#ice vs raw
ggplot()+geom_raster(data=oppar$ice, aes(begin1,begin2,fill=log(N)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+scale_fill_gradient(low="white", high="black")+  theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_ICE_vs_observed.png"), width=10, height=9)+  theme(legend.position = "none") 
#cs vs ice
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(1238*normalized)))+geom_raster(data=oppar$ice, aes(begin2,begin1,fill=log(N)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none") 
#observed vs expected
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(expected)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+
  #geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none") 
#nu
ggplot(melt(data.table(id=biases[,pos],parallel=oppar$par$log_nu),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))
#decay
ggplot(data.table(dist=oppar$binned$distance,decay=oppar$binned$decay)[dist>1e3])+geom_line(aes(dist,decay))+scale_x_log10()+scale_y_log10()
#decay and normalized counts
#oppar$mat[,distance:=pmin(begin2-begin1,4042929+1-(begin2-begin1))]
#setkey(oppar$mat,distance)
#ggplot(data.table(dist=oppar$binned$distance,decay=oppar$binned$decay)[dist>1e3])+geom_line(aes(dist,decay),colour="red")+scale_x_log10()+scale_y_log10()+
#  geom_point(data=oppar$mat,aes(distance,normalized),alpha=0.01)
ggplot(data.table(normalized=counts[,contact.close]*exp(oppar$pred$log_decay_close-oppar$pred$log_mean_cclose),distance=counts[,distance],
                  decay=exp(oppar$pred$log_decay_close))[distance>1000][sample(.N,min(.N,100000))])+
  geom_point(aes(distance,normalized),alpha=0.01)+
  geom_line(aes(distance,decay),colour="red")+scale_x_log10()+scale_y_log10()
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_decay_with_counts.png"), width=10, height=9)
#lfC hist
ggplot(oppar$mat)+geom_histogram(aes(lFC))
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_lFC_hist.png"), width=10, height=9)
#detected interactions and 95% CI for posterior
oppar$mat[,upper:=qgamma(0.975,shape=alpha1,rate=beta1)]
oppar$mat[,lower:=qgamma(0.025,shape=alpha1,rate=beta1)]
oppar$mat[,CI:=upper-lower]
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(normalized)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(CI*normalized/observed)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none")#+scale_colour_brewer(type="div",palette="RdBu")
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_normalized_with_error.png"), width=10, height=9)



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

oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=100000, ncores=30, homogenize=F, ops.count=ops.count, ops.bias=ops.bias)
oppar=postprocess(biases, counts, oppar, resolution=50000, ncores=30, predict.all.means=T)
oppar$ice=iterative_normalization(oppar$mat, niterations=1, resolution=50000, return.binned=T)
save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))


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
