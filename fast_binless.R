library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(scales)

setwd("/Users/yannick/Documents/simulations/cs_norm")


#Welcome to this quick walk through fast binless normalization!
#If you already have a binned matrix with raw counts, go to section "Normalization".
#Otherwise, below is a walk through the whole procedure
#from mapped paired-end reads (in TADbit tsv format) to base-resolution (arrow) plots
#to the binned matrices

### Preprocessing

#generate 3 plots to determine dangling.L, dangling.R and dmin and maxlen
#fread warning can be ignored
#be sure to set the proper read length
#with this Rao dataset, we can start with dangling.L=0 dangling.R=3 maxlen=900 read.len=101 and dmin=1000
a=examine_dataset("/Users/yannick/Dropbox (CRG)/Normalization_Pauli/c133e90d3_chr18_region_cis.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=75)


#Create the CSdata objects that contain each dataset, using the just-determined parameters
#Don't forget to set circularize to the genome length if and only if it is a circular genome
#Be sure to correctly set the condition, replicate and enzyme fields.
#save.data=T is only needed if you plan to take a smaller portion of the data, or if you want to visualize the raw reads.
#You could also filter the tsv file beforehand for that.
#two outputs will be created for a given prefix: prefix_csdata.RData and if save.data==T, prefix_csdata_with_data.RData.
#they contain the respective CSdata objects

dsets=c("c133e90d3","e22e868a9","c133e90d3","e22e868a9")
chrs=c("chr4","chr4","chr18","chr18")

foreach (chr=chrs, name=dsets) %do% {
           cat(chr,name,"\n")
           csd=read_and_prepare(paste0("/Users/yannick/Dropbox (CRG)/Normalization_Pauli/",name,"_",chr,"_region_cis.tsv"),
                                paste0("data/pauli_",name,"_",chr), name, "1",
                                enzyme="MboI", name=paste(chr,name), circularize=-1, dangling.L=c(0),
                                dangling.R=c(3), maxlen=600, read.len=75, dmin=1000, save.data=T)
}


#here we plot the raw reads. We need to load the full csdata object, as only the one without the raw reads is returned.
load("data/pauli_c133e90d3_chr4_csdata_with_data.RData")
data=get_raw_reads(csd@data, csd@biases[,min(pos)], csd@biases[,max(pos)])
plot_binned(data, resolution=10000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])
plot_raw(data[rbegin1>min(rbegin1)+50000&rbegin2<min(rbegin1)+60000])

#we now process the csd objects in a format amenable for the fast binless algorithm
#first, we merge all datasets we want normalized together (they must represent the same genomic region)
#in fast binless, one decay is modelled for all datasets. Thus, all arguments to merge_cs_norm_datasets are ignored.
load("data/pauli_c133e90d3_chr4_csdata.RData")
csd1=csd
load("data/pauli_e22e868a9_chr4_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
#now we bin the raw data at the base resolution we want, and put it in a data table
mat=csnorm:::bin_data(cs,resolution=5000)


### Normalization

#We start with a mat data.table that contains 4 columns. The first three columns are factors
# (or contiguous integers starting at 1) and refer to
# - name: the name of the dataset (remember we normalize multiple datasets together)
# - bin1,bin2: where the counts fall into. We must have bin2 >= bin1
# The last column is observed: how many paired-end reads are observed at these coordinates

#in the fast binless mode, you must ensure no counter diagonal nor row/column is completely zero
#if there are not too many, you can add 1 to the observed counts
rbind(mat[,.(bin=bin1,observed)],mat[,.(bin=bin2,observed)])[,.(sum(observed)),by=bin][V1==0]
mat[,.(d=unclass(bin2)-unclass(bin1),observed)][,sum(observed),by=d][V1==0]
mat[unclass(bin2)-unclass(bin1)>=555,
    observed:=observed+as.integer(c(1,rep(0,length(observed)-1))),
    by=unclass(bin2)-unclass(bin1)]

#the fast binless algorithm computes a binless normalization without estimating
#the fusion penalty and the significance threshold.
#We do a maximum of nouter steps (or less if the relative precision is lower than tol_val)
#and hold the fusion penalty fixed at lam2. Play with lam2 to see its effect.
nouter=20
lam2=10
tol_val=1e-1
bg_steps=5
out=csnorm:::fast_binless(mat, mat[,nlevels(bin1)], lam2, nouter, tol_val, bg_steps)


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



#now we compute differences between the two datasets
ref=1
lam2=10
tol_val=1e-1
diff=as.data.table(csnorm:::fast_binless_difference(out, lam2, ref, tol_val))

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
ggplot(diff)+geom_raster(aes(bin1,bin2,fill=delta))+
  geom_raster(aes(bin2,bin1,fill=delta))+facet_wrap(~name)+
  scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()

