library(binless)
library(data.table)
library(ggplot2)
library(scales)
library(profvis)


#normalize
resolution=20000
fname="data/rao_HiCall_SELP_150k_bpk30_csnorm_optimized.RData"
fname="data/rao_HiCall_FOXP1_1.3M_bpk30_csnorm_optimized.RData"
ncores=30
load(fname)

stuff = binless:::bin_and_chunk(cs, resolution, "all", 1)
counts = binless:::fill_zeros(stuff$counts, stuff$biases, circularize=F, dmin=cs@settings$dmin)
counts = counts[,.(name,id1,id2,pos1,pos2,contact.close,contact.down,contact.far,contact.up,distance)]
counts=merge(counts,stuff$biases[,.(name,id,pos,bin,ibin)],by.x=c("name","id1","pos1"),by.y=c("name","id","pos"))
counts=merge(counts,stuff$biases[,.(name,id,pos,bin,ibin)],by.x=c("name","id2","pos2"),by.y=c("name","id","pos"),suffixes=c("1","2"))
profvis({binless:::estimate_signal(cs,counts,stuff$groups)})

#compare old and new matrix binnings
mat.ori = binless:::predict_binned(cs, resolution, "all", ncores=ncores)
mat.ori[,ori:="old"]
mat.new = binless:::predict_binned_irls(cs, resolution, "all", ncores=ncores)
mat.new[,ori:="new"]
mat=rbind(mat.new,mat.ori)
ggplot(dcast(mat[,.(name,bin1,bin2,ori,log(signal))],name+bin1+bin2~ori))+stat_function(fun=identity)+
  geom_point(aes(new,old,colour=name))
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=signal))+facet_grid(name~ori)
ggplot(mat[,.(name,bin1,bin2,signal,ori,ratio=observed/expected)])+geom_point(aes(signal,ratio))+
  facet_grid(name~ori)+stat_function(fun=identity)

#compare interaction detection
mat.ori = binless:::detect_binned(cs, resolution, "all", ref="expected", threshold=0.95, prior.sd=5, ncores=ncores)
mat.ori[,ori:="old"]
mat.new = binless:::detect_binned_interactions_irls(cs, resolution, "all", threshold=0.95, prior.sd=5, ncores=ncores)
mat.new[,ori:="new"]
mat=rbind(mat.new,mat.ori)
ggplot(dcast(mat[,.(name,bin1,bin2,ori,prob.gt.expected)],name+bin1+bin2~ori))+stat_function(fun=identity)+
  geom_point(aes(new,old,colour=name))
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=log(signal.signif)))+facet_grid(name~ori)+
  geom_point(aes(bin2,bin1,colour=direction),data=mat[is.significant==T])+scale_fill_gradient2()
dcast(mat[,.(name,bin1,bin2,ori,is.significant)],name+bin1+bin2~ori)[,.N,by=c("new","old")]
dcast(mat[,.(name,bin1,bin2,ori,direction)],name+bin1+bin2~ori)[,.N,by=c("new","old")]
mat[is.na(direction)]

#compare difference detection
mat.ori = binless:::detect_binned(cs, resolution, "all", ref="GM MboI 1", threshold=0.95, prior.sd=5, ncores=ncores)
mat.ori[,ori:="old"]
mat.new = binless:::detect_binned_differences_irls(cs, "GM MboI 1", resolution, "all", threshold=0.95, prior.sd=5, ncores=ncores)
mat.new[,ori:="new"]
mat=rbind(mat.new,mat.ori)
ggplot(dcast(mat[,.(name,bin1,bin2,ori,difference)],name+bin1+bin2~ori))+stat_function(fun=identity)+
  geom_point(aes(new,old,colour=name))
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=log(difference.sd)))+facet_grid(name~ori)+
  geom_point(aes(bin2,bin1,colour=direction),data=mat[is.significant==T])+scale_fill_gradient2()
dcast(mat[,.(name,bin1,bin2,ori,is.significant)],name+bin1+bin2~ori)[,.N,by=c("new","old")]
dcast(mat[,.(name,bin1,bin2,ori,direction)],name+bin1+bin2~ori)[,.N,by=c("new","old")]
mat[is.na(direction)]

