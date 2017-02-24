library(csnorm)
library(data.table)
library(ggplot2)
library(scales)
library(profvis)


#normalize
resolution=20000
fname="data/rao_HiCall_SELP_150k_bpk30_csnorm_optimized.RData"
ncores=30
load(fname)

stuff = csnorm:::bin_and_chunk(cs, resolution, "all", 1)
counts = csnorm:::fill_zeros(stuff$counts, stuff$biases, circularize=F, dmin=cs@settings$dmin)
counts = counts[,.(name,id1,id2,pos1,pos2,contact.close,contact.down,contact.far,contact.up,distance)]
counts=merge(counts,stuff$biases[,.(name,id,pos,bin,ibin)],by.x=c("name","id1","pos1"),by.y=c("name","id","pos"))
counts=merge(counts,stuff$biases[,.(name,id,pos,bin,ibin)],by.x=c("name","id2","pos2"),by.y=c("name","id","pos"),suffixes=c("1","2"))
profvis({csnorm:::estimate_signal(cs,counts,stuff$groups)})
