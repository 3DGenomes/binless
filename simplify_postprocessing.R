library(csnorm)
library(data.table)
library(ggplot2)
library(scales)
library(profvis)


#normalize
resolution=5000
fname="data/rao_HiCall_SELP_150k_bpk30_csnorm_optimized.RData"
ncores=30
load(fname)

stuff = csnorm:::bin_and_chunk(cs, resolution, "all", 1)
counts = csnorm:::fill_zeros(stuff$counts, stuff$biases, circularize=F, dmin=cs@settings$dmin)
counts = counts[,.(name,id1,id2,pos1,pos2,contact.close,contact.down,contact.far,contact.up,distance)]
counts=merge(counts,stuff$biases[,.(name,id,pos,bin,ibin)],by.x=c("name","id1","pos1"),by.y=c("name","id","pos"))
counts=merge(counts,stuff$biases[,.(name,id,pos,bin,ibin)],by.x=c("name","id2","pos2"),by.y=c("name","id","pos"),suffixes=c("1","2"))
profvis({csnorm:::estimate_signal(cs,counts,stuff$groups)})

#compare old and new matrix binnings
mat.ori = csnorm:::csnorm_predict_binned(cs, resolution, "all", ncores=ncores)
mat.ori[,ori:="old"]
mat.new = csnorm:::csnorm_predict_binned_irls(cs, resolution, "all", ncores=ncores)
mat.new[,ori:="new"]
mat=rbind(mat.new,mat.ori)
ggplot(dcast(mat[,.(name,bin1,bin2,ori,log(expected.sd))],name+bin1+bin2~ori))+stat_function(fun=identity)+
  geom_point(aes(new,old,colour=name))
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=signal))+facet_grid(name~ori)
ggplot(mat[,.(name,bin1,bin2,signal,ori,ratio=observed/expected)])+geom_point(aes(signal,ratio))+
  facet_grid(name~ori)+stat_function(fun=identity)



