library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


### Preprocessing

#generate 3 plots to determine dangling.L, dangling.R and maxlen respectively (fread warning can be ignored)
#with this Caulobacter dataset, we can start with dangling.L=0 dangling.R=3 maxlen=600 read.len=40 and dmin=2000
a=examine_dataset("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                  skip=0L,nrows=1000000, skip.fbm=F)
#Since the cut site is palindromic, it is expected that the peaks have the same height. You can test this hypothesis.
chisq.test(c(a$data[(rbegin1-re.closest1)==0&strand1==1,.N],
             a$data[(rbegin2-re.closest2)==3&strand2==0,.N]))



#Create the CSdata objects that contain each dataset, using the just-determined parameters
#Don't forget to set circularize to the genome length if and only if it is a circular genome
#be sure to correctly set the condition, replicate and enzyme fields.
#save.data=T is only needed if you plan to take a smaller portion of the data, or if you want to visualize the raw reads.
#You could also filter the tsv file beforehand for that.
#two outputs will be created for a given prefix: prefix_csdata.RData and if save.data==T, prefix_csdata_with_data.RData.
#they contain the respective CSdata objects
csd0=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                      "data/caulo_NcoI_all", "WT", "1", enzyme="NcoI", circularize=4042929, dangling.L=c(0),
                      dangling.R=c(3), maxlen=600, read.len=40, dmin=1000, save.data=T)
csd1=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                      "data/caulo_BglIIr1_all", "WT", "1", enzyme="BglII", circularize=4042929, dangling.L=c(0),
                      dangling.R=c(3), maxlen=600, read.len=40, dmin=1000, save.data=T)
csd2=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate2_reads_int.tsv",
                      "data/caulo_BglIIr2_all", "WT", "2", enzyme="BglII", circularize=4042929, dangling.L=c(0),
                      dangling.R=c(3), maxlen=600, read.len=40, dmin=1000, save.data=T)
csd3=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_rifampicin_reads_int.tsv",
                      "data/caulo_BglIIrif_all", "rif", "1", enzyme="BglII", circularize=4042929, dangling.L=c(0),
                      dangling.R=c(3), maxlen=600, read.len=40, dmin=1000, save.data=T)
csd4=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_novobiocin_reads_int.tsv",
                      "data/caulo_BglIInov_all", "nov", "1", enzyme="BglII", circularize=4042929, dangling.L=c(0),
                      dangling.R=c(3), maxlen=600, read.len=40, dmin=1000, save.data=T)
csd5=read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_dSMC_reads_int.tsv",
                      "data/caulo_BglIIdSMC_all", "dSMC", "1", enzyme="BglII", circularize=4042929, dangling.L=c(0),
                      dangling.R=c(3), maxlen=600, read.len=40, dmin=1000, save.data=T)

#We can plot a particular region of interest
#to save memory, only the CSdata object that does not contain the raw reads is returned.
#We therefore need to load it
load("data/caulo_BglIIr1_all_csdata.RData")
data=get_raw_reads(csd@data, 1080000, 1100000)
plot_binned(data, resolution=1000, b1=1081000, e1=1099000)
plot_raw(data, b1=1081000, e1=1099000)

plot_binned(data, resolution=1000, b1=1089000, e1=1092000)
plot_binned(data, resolution=100, b1=1089000, e1=1092000)
plot_raw(data, b1=1089000, e1=1094000)

### Normalization

#Now it is necessary to create the main CSnorm object. It is a class that holds everything related to this normalization.
#In cut-site normalization, datasets are normalized together, not individually, It is therefore needed to specify how different
#datasets relate to each other. Different enzymes must be modelled by different genomic biases, and this is done automatically 
#if you set the enzyme field properly. Whether all datasets should have the same diagonal decay or not is a more complicated question
#and will vary depending on your experimental design. Here, we want to assess the changes that rifampicin induces, so we use the same
#diagonal decay for both conditions. However, BglII cuts in very different locations than NcoI, and it is expected that the decay be
#a bit different between the two enzymes. As we do not want to detect this, we use one decay for each enzyme.
#Once the cs object is created, have a look at it by typing cs in the interactive prompt, and look for the experimental design matrix.
#The cs object is full of information, be sure to check it out at every step in this script.
load("data/caulo_NcoI_all_csdata.RData")
csd1=csd
load("data/caulo_BglIIr1_all_csdata.RData")
csd2=csd
load("data/caulo_BglIIr2_all_csdata.RData")
csd3=csd
load("data/caulo_BglIIrif_150k_csdata.RData")
csd4=csd
load("data/caulo_BglIInov_all_csdata.RData")
csd5=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3), different.decays="all")
save(cs, file="data/caulo_rif_150k_csnorm.RData")

#The normalization is performed with a fast approximate Gibbs sampler to the full model
#bf_per_kb: number of spline basis functions per genomic kilobase. The larger, the better the fit, but the harder the optimization.
#  as a rule of thumb, start with 0.25 for 6-cutters, and 1 for 4-cutters.
#lambdas: Starting value for genomic spline penalties. Each value of lambda will be run on a different processor, and only the best
#  run will be kept. The value of lambda does not mean anything per se, and it changes with bf_per_kb and the experiment's sequencing
#  depth. Start with 3 to 5 lambda values spaced out between 0.01 and 10.
#ngibbs: number of gibbs sampler iterations to perform. One iteration is often good enough.
#Once optimized, be sure to look at cs@diagnostics. It contains the outputs of the underlying optimization rounds. Make sure that
#in the last runtime.bias and runtime.decay, you see "Convergence detected". Otherwise, the normalization was not successful.
#You might either want to increase ngibbs or iter, or both.
load("data/caulo_rif_150k_csnorm.RData")
cs = normalize_binless(cs, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, ngibbs = 1, iter=10000, init_alpha=1e-7, ncounts = 1000000)
save(cs, file="data/caulo_rif_csnorm_optimized.RData")

#Some diagnostics plots can be generated. The function generates 4 plots and also returns a data table.
#The p-values of most points reported should not be too extreme. If they are, then the fit went bad.
check=check_fit(cs)
#in particular, the following should be close to 0.05
check$counts[pval<0.05,.N]/check$counts[,.N]
check$counts[pval>0.95,.N]/check$counts[,.N]

nu=data.table(name=cs@biases[,name],pos=cs@biases[,pos],nu=exp(cs@par$log_nu),delta=exp(cs@par$log_delta))
ggplot(nu)+geom_line(aes(pos,log(nu)))+facet_wrap(~name,scales = "free")

ggplot(cs@par$biases)+geom_line(aes(pos,log_nu))+
  geom_pointrange(aes(pos,log_nu+z,ymin=log_nu+z+std,ymax=log_nu+z-std,colour=cat), alpha=0.1)+
  facet_wrap(~name)+ylim(-5,5)

ggplot(cs@par$decay[dist>1000])+geom_line(aes(dist,log_decay),colour="red")+ facet_wrap(~name) +
  geom_pointrange(aes(x=dist,y=log_decay+z,ymin=log_decay+z-std,ymax=log_decay+z+std),alpha=0.1)+
  scale_x_log10()

ggplot(dmat[dist>1000])+geom_line(aes(dist,log_decay),colour="red")+ facet_wrap(~name) +
  geom_pointrange(aes(x=dist,y=log_decay+z,ymin=log_decay+z-std,ymax=log_decay+z+std),alpha=0.1)+
  scale_x_log10()


### Binning at a given resolution

#Once normalized, we generate binned matrices to perform the interaction detection. Here we choose to bin at 20kb resolution.
#there are three dispersion calculations in this beta version, you might want to try them all out.
#You can optionally ask for computation of ICE-normalized matrices, for comparison.
load("data/caulo_rif_csnorm_optimized.RData")
cs=bin_all_datasets(cs, resolution=20000, ncores=30, verbose=T, ice=1)
save(cs, file="data/caulo_rif_csnorm_optimized.RData")

#To get the binned matrices, use the following command with the same arguments than those passed to bin_all_datasets
#Since you did not do any grouping yet (see below), pass group="all".
mat=get_matrices(cs, resolution=20000, group="all")
#The mat data.table contains a number of matrices:
#observed: the raw counts in each bin
#expected and expected.sd: the expected number of counts, according to csnorm, along with their standard deviation
#normalized and normalized.sd: The matrix of observed/expected, with standard deviations
#decaymat: the diagonal decay matrix
#ice.x: If computed, an ICE matrix with x iterations
#icelike and icelike.sd: A "normalized" matrix with diagonal decay, useful to compare with ICEd matrices, among others

#for example, we can plot the ICE-like matrix in the upper corner, with corresponding error bars in the lower corner
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike.sd)))+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

#For that binning, we can also verify whether the normalization has been performed correctly, using a residual plot
ggplot(mat)+geom_point(aes(expected,observed),alpha=0.1)+stat_function(fun=identity,colour="red")+facet_wrap(~name)+
  scale_x_log10(limits=c(0.5,NA))+scale_y_log10(limits=c(0.5,NA))
#or a density plot
ggplot(mat)+geom_density(aes(log10(normalized),colour=name))+facet_wrap(~name)




### Interaction calling 

#To detect interactions, you need to specify an already binned dataset by providing resolution and detection.type.
#Since you did not do any grouping yet (see below), pass group="all".
#The interaction detection is made at a 95% posterior confidence threshold
cs=detect_interactions(cs, resolution=20000, group="all", threshold=0.95, ncores=30)
save(cs,file="data/caulo_BglII_decay_all_csnorm_optimized.RData")

#you can view the called interactions
#since there was no grouping, pass group="all"
#for simple interactions pass type="interactions" and ref="expected"
mat=get_interactions(cs, type="interactions", resolution=20000, group="all",
                     threshold=0.95, ref="expected")
#for example, we can plot the ice-like matrices with highlighted interactions in the upper left corner
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_grid(~name)

cs=detect_differences(cs, ref="WT BglII 2", resolution=20000, group="all", threshold=0.95, ncores=30)
dmat=get_interactions(cs, type="differences", resolution=20000, group="all",
                     threshold=0.95, ref="WT BglII 2")
ggplot(dmat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_grid(~name)


### Grouping of datasets

#Datasets can be grouped (merged) to improve the quality of the matrices and the interaction detection power
#for example, we can group datasets by condition and enzyme. For that purpose, pass group the corresponding groupings.
cs=group_datasets(cs, resolution=20000, group=c("condition", "enzyme"), ice=1, verbose=T)

#Again, you can detect interactions in these grouped datasets, and plot the results
cs=detect_interactions(cs, resolution=20000, group=c("condition", "enzyme"), detection.type=1,
                        threshold=0.95, ncores=30)
mat=get_interactions(cs, type="interactions", resolution=20000, group=c("condition", "enzyme"), detection.type=1,
                     threshold=0.95, ref="expected")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)


### Difference calling

#you can call significant differences. It's just like calling interactions, but you need to specify a reference.
#All other matrices will then be compared to that one.
cs=detect_differences(cs, ref="WT BglII 2", resolution=20000, group="all", threshold=0.95, ncores=30)
save(cs, file="data/caulo_BglII_decay_all_csnorm_optimized.RData")

#the interactions can be retrieved as usual, with a few changes
# specify type="differences" and ref as given in detect_differences
mat=get_interactions(cs, type="differences", resolution=20000, group="all",
                     threshold=0.95, ref="WT BglII 2")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient(low="white", high="black", na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

#You can also plot the matrix of the ratio of one experiment to its reference, along with its error bar
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=log(ratio)))+
  geom_raster(aes(begin2,begin1,fill=log(ratio.sd)))+
  scale_fill_gradient2(na.value = "white")+theme_bw()+theme(legend.position = "none")+
  facet_wrap(~name)

save(cs, file="data/caulo_csnorm_optimized.RData")



### END


#nu and delta correlation
nu = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  binless:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,dset=j)]
}
nu[,pbin:=cut(pos,3)]
ggplot(nu)+geom_line(aes(pos,exp(log_nu),colour=dset))+facet_wrap(~pbin,scales = "free_x", nrow=3)+
  #scale_y_continuous(limits = c(0,2))+
  ylab("nu")+scale_y_log10()
#ggsave(filename=paste0("images/",prefix,"_nu.png"), width=10, height=7.5)

#sorted with data points
nu.obs = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  a=data.table(id=cs@counts[,id1], pos=cs@counts[,pos1],
               obs=cs@counts[,contact.down]*exp(-cs@pred$log_mean_cdown),
               dset=j)[sample(.N,min(.N,100000))]
  a=merge(a,cs@biases[,.(id,nu=exp(cs@par$log_nu))],by="id")
  a[,obs:=obs*nu]
  setkey(a,nu,pos)
  a[,rank:=.I]
  a
}
ggplot(nu.obs)+geom_line(aes(rank,nu,colour=dset))+geom_point(aes(rank,obs,colour=dset),alpha=0.01)+
  #scale_y_continuous(limits = c(0,2))+
  ylab("nu")+scale_y_log10()


#diagonal decay
decay.obs = foreach (i=fnames,j=dsets,.combine=rbind) %do% {
  load(i)
  data.table(dist=cs@counts[,distance],
             obs=cs@counts[,contact.down]*exp(-cs@pred$log_mean_cdown+cs@pred$log_decay_down),
             decay=exp(cs@pred$log_decay_down),
             dset=j)[sample(.N,min(.N,100000))]
}
decay.obs=decay.obs[dist>1000]
setkey(decay.obs, dist)
ggplot(decay.obs)+geom_line(aes(dist,decay,colour=dset))+geom_point(aes(dist,obs,colour=dset),alpha=0.01)+
  scale_x_log10()+scale_y_log10()
#ggsave(filename=paste0("images/",prefix,"_decay.png"), width=10, height=7.5)

