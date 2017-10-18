library(ggplot2)
library(scales)
library(data.table)
library(binless)


#Welcome to this quick walk through fast binless normalization!
#If you already have a binned matrix with raw counts, go to section "Normalization".
#Otherwise, below is a walk through the whole procedure
#from mapped paired-end reads (in TADbit tsv format) to base-resolution (arrow) plots
#to the binned matrices

### Preprocessing

#The following command returns a list that contains three plots and the data used to generate them
#These plots are useful to determine dangling.L, dangling.R, dmin and maxlen
#filename is passed to fread, so you can use zcat to read compressed files directly
#fread warning can be ignored
#be sure to set the proper read length
#with this Rao dataset, we can start with dangling.L=0 dangling.R=3 maxlen=900 read.len=101 and dmin=1000
a=examine_dataset("zcat example/GM12878_MboI_HICall_FOXP1ext.tsv.gz",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)


#Create the CSdata objects that contain each dataset, using the just-determined parameters
#Fast binless does not support circular genomes, so leave circular=-1
#Be sure to correctly set the condition, replicate and enzyme fields.
#save.data=T is only needed if you plan to take a smaller portion of the data, or if you want to visualize the raw reads.
#You could also filter the tsv file beforehand for that.
#two outputs will be created for a given prefix: prefix_csdata.RData and if save.data==T, prefix_csdata_with_data.RData.
#they contain the respective CSdata objects.

name="FOXP1ext"
size="2.3M"
dsets=c("GM12878","IMR90")
maxlens=c(900,600)
foreach (dset=dsets, maxlen=maxlens) %do% {
  csd=read_and_prepare(paste0("zcat example/",dset,"_MboI_HICall_",name,".tsv.gz"),
                       paste0("example/rao_HiCall_",dset,"_",name,"_",size), dset, "1",
                       enzyme="MboI", name=paste(name,dset,"all"), circularize=-1, dangling.L=c(0),
                       dangling.R=c(3), maxlen=maxlen, read.len=101, dmin=1000, save.data=T)
}         


#here we plot the raw (unnormalized) reads. We need to load the full csdata object,
#as only the one without the raw data is returned.

load("example/rao_HiCall_GM12878_FOXP1ext_2.3M_csdata_with_data.RData")
data=get_raw_reads(csd@data, csd@biases[,min(pos)], csd@biases[,max(pos)])
#plot the whole 2.3M region at 10kb
plot_binned(data, resolution=10000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])
#plot a 20kb subset of it with base resolution (arrow plot)
plot_raw(data, b1=data[,min(rbegin1)+50000], e1=data[,min(rbegin1)+70000])


#We now process the csd objects in a format amenable for the fast binless algorithm
#first, we merge all datasets we want normalized together (they must represent the same genomic region)
#in fast binless, one decay is modelled for all datasets. Thus, all arguments to merge_cs_norm_datasets are ignored.
load("example/rao_HiCall_GM12878_FOXP1ext_2.3M_csdata.RData")
csd1=csd
load("example/rao_HiCall_IMR90_FOXP1ext_2.3M_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
#now we bin the raw data at the base resolution we want, and put it in a data table
mat=binless:::bin_data(cs,resolution=5000)
#write.table(mat,file = "example/rao_HiCall_FOXP1ext_2.3M_mat_5kb.dat", quote = T, row.names = F)

### Normalization

#mat=fread("zcat example/rao_HiCall_FOXP1ext_2.3M_mat_5kb.dat.gz",stringsAsFactors = T)

#We start with a mat data.table that contains 4 columns. The first three columns are factors
# (or contiguous integers starting at 1) and refer to
# - name: the name of the dataset (remember we normalize multiple datasets together)
# - bin1,bin2: where the counts fall into. We must have bin2 >= bin1
# The last column is observed: how many paired-end reads are observed at these coordinates

#in the fast binless mode, you must ensure no counter diagonal nor row/column is completely zero
#if there are not too many, you can add 1 to the observed counts
#also, because of the simplified decay, the very last bin cannot be zero in all datasets
rbind(mat[,.(bin=bin1,observed)],mat[,.(bin=bin2,observed)])[,.(sum(observed)),by=bin][V1==0]
mat[,.(d=unclass(bin2)-unclass(bin1),observed)][,sum(observed),by=d][V1==0]
mat[unclass(bin1)==1&unclass(bin2)==max(unclass(bin2))]

#the fast binless algorithm computes a binless normalization without estimating
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
