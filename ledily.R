library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(scales)

setwd("/home/yannick/simulations/cs_norm")


### Preprocessing

a=examine_dataset("/scratch/ledily/mapped/T47D_estrogen_60/capture_chr6_149579000-153679800.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=125)
a=examine_dataset("/scratch/ledily/mapped/T47D_Progesterone_30/capture_chr6_149579000-153679800.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=125)
a=examine_dataset("/scratch/ledily/mapped/T47D_Progesterone_180/capture_chr6_149579000-153679800.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=125)

fnames=c("/scratch/ledily/mapped/T47D_estrogen_60/capture_chr6_149579000-153679800.tsv",
         "/scratch/ledily/mapped/T47D_Progesterone_30/capture_chr6_149579000-153679800.tsv",
         "/scratch/ledily/mapped/T47D_Progesterone_180/capture_chr6_149579000-153679800.tsv")
prefixes=c("data/ledily_T47D_es_60", "data/ledily_T47D_pg_30", "data/ledily_T47D_pg_180")
names=c("T47D es 60", "T47D pg 30", "T47D pg 180")

registerDoParallel(cores=30)
foreach (fname=fnames,prefix=prefixes,name=names) %dopar%
  read_and_prepare(fname, prefix, name, "1",
                   enzyme="MboI", circularize=-1, dangling.L=c(0),
                   dangling.R=c(3), maxlen=750, read.len=125, dmin=1000, save.data=T)

 
#plot raw reads
load("data/ledily_T47D_es_60_csdata_with_data.RData")
data=get_raw_reads(csd@data, 150065065-10000, 150065065+10000)
plot_binned(data, resolution=1000, b1=150065065-10000, e1=150065065+10000)
plot_raw(data, b1=150065065-10000, e1=150065065+10000)

#plot approximate coverage
load("data/ledily_T47D_es_60_csdata.RData")
ggplot(csd@biases)+geom_point(aes(pos,rejoined))
ggplot(csd@biases)+geom_point(aes(pos,rejoined))+xlim(149570000,149590000)
ggplot(csd@biases)+geom_point(aes(pos,rejoined))+xlim(153675000,153685000)

#determine new boundaries
csd@biases[pos>149570000&pos<149590000] #lower boundary is 149579000
csd@biases[pos>153675000&pos<153685000] #upper boundary is 153679800

load("data/ledily_T47D_estrogen_60_csdata.RData")
csd1=csd
load("data/ledily_T47D_progesterone_30_csdata.RData")
csd2=csd
cs=merge_cs_norm_datasets(list(csd1,csd2), different.decays = "none")
cs = run_gauss(cs, bf_per_kb=3, bf_per_decade=10, bins_per_bf=10, ngibbs = 10, iter=100000, init_alpha=1e-7, ncounts = 1000000)
save(cs, file="data/ledily_T47D_E60P30_csnorm_optimized.RData")

#Once optimized, be sure to look at cs@diagnostics$plot. It monitors convergence of the normalization. Make sure that
#all three panels reach a plateau, and the last points report "Convergence detected". Otherwise, the normalization was not successful.
#You might then either want to increase ngibbs or iter, or both.
cs@diagnostics$plot
#similar plots can be generated. For example, to visualize alpha
ggplot(cs@diagnostics$params[leg=="disp",.(alpha=sapply(alpha,function(x){x[[1]]}),out.last),by=step])+
  geom_line(aes(step,alpha))+geom_point(aes(step,alpha,colour=out.last))
#or eC
ggplot(cs@diagnostics$params[,.(step,leg,eC=sapply(eC,function(x){x[[1]]}),out.last)])+
  geom_line(aes(step,eC))+geom_point(aes(step,eC,colour=out.last))+facet_grid(~leg)
ggplot(cs@diagnostics$params[,.(step=.I/3,leg,eC=sapply(eC,function(x){x[[1]]}),out.last)])+
  geom_line(aes(step,eC))+geom_point(aes(step,eC,colour=leg))

#To assess normalization quality, some other diagnostics plots can be generated as well.
#The function generates 4 plots and also returns a data table.
#The p-values of most points reported should not be too extreme. If they are, then the fit went bad.
check=check_fit(cs)
#in particular, the following should be close to 0.05
check$counts[pval<0.05,.N]/check$counts[,.N]
check$counts[pval>0.95,.N]/check$counts[,.N]


#you can also plot the value of the biases along the genome. Here, we plot nu and the reduced aggregate counts used for optimization.
ggplot(cs@par$biases)+geom_pointrange(aes(pos,etahat,ymin=etahat-std,ymax=etahat+std,colour=cat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_grid(name ~ cat)#+
xlim(73800000,73810000)+ylim(-10,10)
xlim(26627203,26637203)+ylim(-10,10)


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
  geom_raster(aes(begin1,begin2,fill=log(icelike)))+
  geom_raster(aes(begin2,begin1,fill=log(icelike)))+
  scale_fill_gradient(low="white", high="black")+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  geom_raster(aes(begin2,begin1,fill=-log(normalized)))+
  scale_fill_gradient2()+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name)
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=observed))+
  geom_raster(aes(begin2,begin1,fill=expected))+
  scale_fill_gradient(low="white", high="black")+
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


cs=detect_differences(cs, resolution=5000, group="all", threshold=0.95, ncores=30, ref="GM MboI 1")
mat=get_interactions(cs, type="differences", resolution=5000, group="all", threshold=0.95, ref="GM MboI 1")
ggplot(mat)+
  geom_raster(aes(begin1,begin2,fill=-log(normalized)))+
  #geom_raster(aes(begin2,begin1,fill=-log(signal)))+
  geom_point(aes(begin2,begin1,colour=direction),data=mat[is.significant==T])+
  scale_fill_gradient2()+ scale_colour_manual(values = muted(c("blue","red")))+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~name, nrow=2, ncol=2)



