library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(scales)

setwd("/home/yannick/simulations/cs_norm")


### Preprocessing

#generate 3 plots to determine dangling.L, dangling.R and maxlen respectively (fread warning can be ignored)
#with this Rao dataset, we can start with dangling.L=0 dangling.R=3 maxlen=900 read.len=101 and dmin=1000
a=examine_dataset("/scratch/rao/mapped/GM12878_MboI_in_situ/HICall_both_filled_map_chrX_Peak1.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
a=examine_dataset("/scratch/rao/mapped/GM12878_MboI_in_situ/HICall_both_filled_map_chr21_ADAMT.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
a=examine_dataset("/scratch/rao/mapped/IMR90_MboI_in_situ/HICall_IMR90_both_map_chr21_ADAMT.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)


#Create the CSdata objects that contain each dataset, using the just-determined parameters
#Don't forget to set circularize to the genome length if and only if it is a circular genome
#Be sure to correctly set the condition, replicate and enzyme fields.
#save.data=T is only needed if you plan to take a smaller portion of the data, or if you want to visualize the raw reads.
#You could also filter the tsv file beforehand for that.
#two outputs will be created for a given prefix: prefix_csdata.RData and if save.data==T, prefix_csdata_with_data.RData.
#they contain the respective CSdata objects

foreach (chr=c("chrX","chr1","chr21","chr1","chr8","chr3","chr4","chr21"),
         name=c("Peak1","SELP","ADAMT","Talk","SEMA3C","FOXP1","PARM1","Comparison"),
         size=c("450k","150k","2.3M","2M","1M","1.3M","600k","1.7M")) %do% {
           csd=read_and_prepare(paste0("/scratch/rao/mapped/GM12878_MboI_in_situ/HICall_both_filled_map_",chr,"_",name,".tsv"),
                                paste0("data/rao_HiCall_GM12878_",name,"_",size), "GM", "1",
                                enzyme="MboI", circularize=-1, dangling.L=c(0),
                                dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)
           csd=read_and_prepare(paste0("/scratch/rao/mapped/IMR90_MboI_in_situ/HICall_IMR90_both_map_",chr,"_",name,".tsv"),
                                paste0("data/rao_HiCall_IMR90_",name,"_",size), "IMR90", "1", enzyme="MboI", circularize=-1, dangling.L=c(0),
                                dangling.R=c(3), maxlen=600, read.len=101, dmin=1000, save.data=T)
         }


#here we plot the raw reads. We need to load the full csdata object, as only the one without the raw reads is returned.
load("data/rao_HiCall_chrX_450k_csdata_with_data.RData")
data=get_raw_reads(csd@data, csd@biases[,min(pos)], csd@biases[,max(pos)])
plot_binned(data, resolution=10000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])


#The normalization is performed with a fast approximation to the full model
#bf_per_kb: number of spline basis functions per genomic kilobase. The larger, the better the fit, but the harder the optimization.
#  as a rule of thumb, start with 1 for 6-cutters, and 5 for 4-cutters.
#lambdas: Starting value for genomic spline penalties. Each value of lambda will be run on a different processor, and only the best
#  run will be kept. The value of lambda does not mean anything per se, and it changes with bf_per_kb and the experiment's sequencing
#  depth. Start with 10 lambda values spaced out between 0.01 and 100.
#ngibbs: number of gibbs sampler iterations to perform. Try 20 for a start.
#prefix: if you set this, individual optimizations will be stored using that prefix. Can be useful for recovery of a failed normalization.
load("data/rao_HiCall_chrX_450k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays = "none")
cs = run_gauss(cs, bf_per_kb=5, bf_per_decade=10, bins_per_bf=10, lambdas=10**seq(from=-2,to=2,length.out=10),
               ngibbs = 20, iter=10000, verbose=T,
               prefix="tmp/rao_HiCall_chrX_450k", ncores=10)
save(cs, file="data/rao_HiCall_chrX_450k_csnorm_optimized.RData")

#Once optimized, be sure to look at cs@diagnostics$plot. It monitors convergence of the normalization. Make sure that
#all three panels reach a plateau, and the last points report "Convergence detected". Otherwise, the normalization was not successful.
#You might then either want to increase ngibbs or iter, or both.
cs@diagnostics$plot
#similar plots can be generated. For example, to visualize lambda_nu
ggplot(cs@diagnostics$params[leg=="bias",.(lambda_nu=sapply(lambda_nu,function(x){x[[1]]}),out.last),by=step])+
  geom_line(aes(step,lambda_nu))+geom_point(aes(step,lambda_nu,colour=out.last))
#or eC
ggplot(cs@diagnostics$params[,.(step,leg,eC=sapply(eC,function(x){x[[1]]}),out.last)])+
  geom_line(aes(step,eC))+geom_point(aes(step,eC,colour=out.last))+facet_grid(~leg)

#To assess normalization quality, some other diagnostics plots can be generated as well.
#The function generates 4 plots and also returns a data table.
#The p-values of most points reported should not be too extreme. If they are, then the fit went bad.
check=check_fit(cs)
#in particular, the following should be close to 0.05
check$counts[pval<0.05,.N]/check$counts[,.N]
check$counts[pval>0.95,.N]/check$counts[,.N]


#you can also plot the value of the biases along the genome. Here, we plot nu and the reduced aggregate counts used for optimization.
ggplot(cs@par$biases)+geom_line(aes(pos,log_nu))+
  geom_pointrange(aes(pos,log_nu+z,ymin=log_nu+z-std,ymax=log_nu+z+std,colour=cat), alpha=0.1)+
  facet_wrap(~name,scales = "free")

#the diagonal decay can be plotted this way
ggplot(cs@par$decay)+geom_line(aes(dist,log_decay))+
  geom_pointrange(aes(dist,log_decay+z,ymin=log_decay+z-std,ymax=log_decay+z+std), alpha=0.1)+
  facet_wrap(~name,scales = "free")+scale_x_log10()




### Binning at a given resolution

#Once normalized, we generate binned matrices to perform the interaction detection. Here we choose to bin at 10kb resolution.
#You can optionally ask for computation of ICE-normalized matrices, for comparison.
load("data/rao_HiCall_chrX_450k_csnorm_optimized.RData")
cs=bin_all_datasets(cs, resolution=5000, ncores=30, verbose=T, ice=100)
save(cs, file="data/rao_HiCall_SELP_150k_csnorm_optimized_newton.RData")

#To get the binned matrices, use the following command with the same arguments than those passed to bin_all_datasets
#Since you will not do any grouping in this example, pass group="all".
mat=get_matrices(cs, resolution=5000, group="all")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=-log(normalized)))+
  scale_fill_gradient2()+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)


#For that binning, we can also verify whether the normalization has been performed correctly, using a residual plot
ggplot(mat)+geom_point(aes(expected,observed),alpha=0.1)+stat_function(fun=identity,colour="red")+facet_wrap(~name)+
  scale_x_log10(limits=c(0.5,NA))+scale_y_log10(limits=c(0.5,NA))
#or a density plot
ggplot(mat)+geom_density(aes(log10(normalized),colour=name))+facet_wrap(~name)

load("data/rao_HiCall_GM12878_bpk1_fix_Peak1_450k_csnorm_optimized.RData")
load("data/rao_HiCall_GM12878_bpk1_variable_Peak1_450k_csnorm_optimized.RData")
load("data/rao_HiCall_GM12878_bpk1_newbiases_variable_Peak1_450k_csnorm_optimized.RData")
cs@diagnostics$plot
ggsave(filename="diag_variable_new.pdf",width=10,height=5)
cs@diagnostics$runtime
cs@par$value


### Interaction calling 

#To detect interactions, you need to specify an already binned dataset by providing resolution and detection.type.
#Since you did not do any grouping, pass group="all".
#The interaction detection is made at a 95% posterior confidence threshold
load("data/rao_HiCall_chrX_450k_csnorm_optimized.RData")
cs=detect_interactions(cs, resolution=5000, group="all", threshold=0.95, ncores=30)
save(cs, file="data/rao_HiCall_chrX_450k_csnorm_optimized.RData")

#you can view the called interactions
#since there was no grouping, pass group="all"
#for simple interactions pass type="interactions" and ref="expected"
mat=get_interactions(cs, type="interactions", resolution=5000, group="all", threshold=0.95, ref="expected")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  #geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name, nrow=2, ncol=2)


