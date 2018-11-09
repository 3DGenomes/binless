library(ggplot2)
library(scales)
library(data.table)
library(foreach)
library(binless)

### Welcome to this quick guide to binless normalization!

### Preprocessing

#Here, we will prepare the data and make it amenable to normalization with binless.
#We will plot the data both in raw (binned) format, and using base-resolution (arrow) plots

#The following command returns a list that contains three plots and the data used to generate them
#These plots are useful to determine dangling.L, dangling.R, dmin and maxlen
#check out the FAQ if there's not enough information in the comments
# https://github.com/3DGenomes/binless/issues?q=label%3AFAQ
#filename is passed to fread, so you can use zcat to read compressed files directly
#fread warning can be ignored
#be sure to set the proper read length
#with this Rao dataset, we can start with dangling.L=0 dangling.R=3 maxlen=900 read.len=101 and dmin=1000
#if you use a Mac, use gzcat instead of zcat, or provide the path to an uncompressed file
#refer to README.md for a description of the tsv file format
a=examine_dataset("zcat example/GM12878_MboI_HICall_SEMA3C.tsv.gz",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
a$data #the data that was read.
#be sure to set the proper read.len, which should be
a$data[,median(length1)]
a$pdangling #select begin of dangling ends. With DpnII, one can expect 0 for dangling.L and 3 for dangling.R
a$pdiag #the distribution of sonication fragments. In this experiment, reads are mostly smaller than maxlen=900
a$pclose #You should see an elbow at short distance, which should not be included. If not, try increasing nrows.
         #Set dmin after it, e.g. dmin=1000. Do not set dmin<maxlen.

#Create the CSdata objects that contain each dataset, using the just-determined parameters
#Fast binless does not support circular genomes, so leave circularize=-1
#Be sure to correctly set the condition, replicate and enzyme fields.
#save.data=T is only needed if you plan to take a smaller portion of the data, or if you want to visualize the raw reads.
#You could also filter the tsv file beforehand for that.
#two outputs will be created for a given prefix: prefix_csdata.RData and if save.data==T, prefix_csdata_with_data.RData.
#they contain the respective CSdata objects.

name="SEMA3C"
size="1M"
dsets=c("GM12878","IMR90")
maxlens=c(900,600)
foreach (dset=dsets, maxlen=maxlens) %do% {
  csd=read_and_prepare(paste0("zcat example/",dset,"_MboI_HICall_",name,".tsv.gz"),
                       paste0("example/rao_HiCall_",dset,"_",name,"_",size), dset, "1",
                       enzyme="MboI", name=paste(name,dset,"all"), circularize=-1, dangling.L=c(0),
                       dangling.R=c(3), maxlen=maxlen, read.len=101, dmin=1000, save.data=T)
}         


#In the previous step, you generated two files, one ending in "_with_data.RData"
#which contains a csdata object that includes the raw reads. Here, we plot them
#using the arrow plot. When you load the file, it creates a csd object that
#contains a csd@data data.table with the necessary info.
load("example/rao_HiCall_GM12878_SEMA3C_1M_csdata_with_data.RData")
data=get_raw_reads(csd@data, csd@biases[,min(pos)], csd@biases[,max(pos)])
#plot the whole 1M region at 5kb resolution
plot_binned(data, resolution=5000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])
#plot a 20kb subset of it with base resolution (arrow plot)
plot_raw(data, b1=data[,min(rbegin1)+50000], e1=data[,min(rbegin1)+70000])



#We now process the csd objects in a format amenable for binless
#we merge all datasets we want normalized together (they must represent the same genomic region)
csd1=get(load("example/rao_HiCall_GM12878_SEMA3C_1M_csdata.RData"))
csd2=get(load("example/rao_HiCall_IMR90_SEMA3C_1M_csdata.RData"))
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays="none")
#save(cs,file="example/rao_HiCall_SEMA3C_csnorm.RData")


