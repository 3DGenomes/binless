library(ggplot2)
library(scales)
library(data.table)
library(binless)

### Normalization using the fast and approximate binless algorithm

#Fast binless requires a binned matrix of raw data at high resolution.
#You can either load it yourself using data.table::fread (file format description in README.md)
#don't forget to use stringsAsFactors=T
#On a mac, use gzcat instead of zcat (or feed it the uncompressed file)
#mat=fread("zcat example/rao_HiCall_FOXP1ext_2.3M_mat_5kb.dat.gz",stringsAsFactors = T)

#You can also extract it from the CSnorm object built during preprocessing
#load("example/rao_HiCall_FOXP1ext_csnorm.RData")
mat=binless:::bin_data(cs,resolution=5000)
#write.table(mat,file = "example/rao_HiCall_FOXP1ext_2.3M_mat_5kb.dat", quote = T, row.names = F)


#In fast binless, you must ensure no counter diagonal nor row/column is completely zero
#if there are not too many, you can add 1 to the observed counts
#also, because of the simplified decay, the very last bin cannot be zero in all datasets
#empty rows:
rbind(mat[,.(bin=bin1,observed)],mat[,.(bin=bin2,observed)])[,.(sum(observed)),by=bin][V1==0]
#empty counter diagonals:
mat[,.(d=unclass(bin2)-unclass(bin1),observed)][,sum(observed),by=d][V1==0]
#last bin:
mat[unclass(bin1)==1&unclass(bin2)==max(unclass(bin2))]

#The fast binless algorithm computes a binless normalization without estimating
#the fusion penalty and the significance threshold.
#We do a maximum of nouter steps (or less if the relative precision is lower than tol_val)
#and hold the fusion penalty fixed at lam2. Play with lam2 to see its effect.
nouter=20
lam2=5
tol_val=2e-1
bg_steps=5
out=binless:::fast_binless(mat, mat[,nlevels(bin1)], lam2, nouter, tol_val, bg_steps)


#Here follow pretty much all the plots you could think of (be sure to check out signal and binless)
#all data
a=as.data.table(out$mat)
#observed matrix
plot_binless_matrix(out$mat, upper="log(observed)", lower="log(observed)")
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=log(observed)))+geom_raster(aes(bin2,bin1,fill=log(observed)))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"), na.value = "white")
#expected matrix
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=log(expected)))+geom_raster(aes(bin2,bin1,fill=log(expected)))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#fitted biases
ggplot(data.table(bin=1:nlevels(mat[,bin1]),log_biases=out$log_biases))+geom_point(aes(bin,log_biases,colour="cpp"))#+
#biases matrix
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=biases))+geom_raster(aes(bin2,bin1,fill=biases))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#fitted decay
ggplot(data.table(distance=1:nlevels(mat[,bin1]),log_decay=out$log_decay))+geom_point(aes(distance,log_decay,colour="cpp"))#+
#decay matrix
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=decay))+geom_raster(aes(bin2,bin1,fill=decay))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#signal matrix
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=signal))+geom_raster(aes(bin2,bin1,fill=signal))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#signal and phihat
ggplot(as.data.table(out$mat))+geom_raster(aes(bin1,bin2,fill=signal))+geom_raster(aes(bin2,bin1,fill=pmin(phihat,max(signal))))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#signal and log(observed)-background
ggplot(as.data.table(out$mat))+geom_raster(aes(bin1,bin2,fill=signal))+geom_raster(aes(bin2,bin1,fill=log(observed/expected)+signal))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"), na.value = "white")
#weights
ggplot(as.data.table(out$mat))+geom_raster(aes(bin1,bin2,fill=log(weights)))+geom_raster(aes(bin2,bin1,fill=log(weights)))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#binless matrix
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=binless))+geom_raster(aes(bin2,bin1,fill=binless))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
#binless and observed
ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=binless))+geom_raster(aes(bin2,bin1,fill=log(observed)*max(binless)/log(max(observed))))+
  facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"), na.value = "white")
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_binless_signal.pdf", width=18,height=8)


#now we compute differences between the two datasets
ref=1
lam2=10
tol_val=2e-1
diff=as.data.table(binless:::fast_binless_difference(out, lam2, ref, tol_val))

#some plots (check out delta)
#differences: log(observed)
ggplot(diff)+geom_raster(aes(bin1,bin2,fill=log(observed)))+
  geom_raster(aes(bin2,bin1,fill=log(observed)))+facet_wrap(~name)+
  scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"),na.value = "white")+coord_fixed()
#differences: phi_ref vs phihat
ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(phi_ref,5)))+
  geom_raster(aes(bin2,bin1,fill=pmin(phihat,5)),data=out$mat)+facet_wrap(~name)+
  scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
#differences: delta
ggplot(diff[name!=ref])+geom_raster(aes(bin1,bin2,fill=delta))+
  geom_raster(aes(bin2,bin1,fill=delta))+facet_wrap(~name)+
  scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_binless_difference.pdf", width=10,height=8)
