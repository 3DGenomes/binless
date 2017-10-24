library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)
library(scales)

setwd("/Users/yannick/Documents/simulations/cs_norm")


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
mat=binless:::bin_data(cs,resolution=2500)

#in the fast binless mode, you must ensure no counter diagonal nor row/column is completely zero
#if there are not too many, you can add 1 to the observed counts
rbind(mat[,.(bin=bin1,observed)],mat[,.(bin=bin2,observed)])[,.(sum(observed)),by=bin][V1==0]
mat[,.(d=unclass(bin2)-unclass(bin1),observed)][,sum(observed),by=d][V1==0]
mat[unclass(bin2)-unclass(bin1)>=513,
    observed:=observed+as.integer(c(1,rep(0,length(observed)-1))),
    by=unclass(bin2)-unclass(bin1)]

#the fast binless algorithm computes a binless normalization without estimating
#the fusion penalty and the significance threshold.
#We do a maximum of nouter steps (or less if the relative precision is lower than tol_val)
#and hold the fusion penalty fixed at lam2. Play with lam2 to see its effect.
nouter=20
lam2=5
tol_val=1e-1
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



#now we compute differences between the two datasets
ref=1
lam2=10
tol_val=1e-1
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
ggplot(diff)+geom_raster(aes(bin1,bin2,fill=delta))+
  geom_raster(aes(bin2,bin1,fill=delta))+facet_wrap(~name)+
  scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()



#cross-validation diagnostics
lam2s=10^seq(-1,1,length.out=20)
#lam2s=10**seq(0.9,-0.2,length.out=10)
#lam2s=c(0,10**seq(0,-10,length.out=11))
ret=binless:::fast_binless_eval_cv(out, lam2s, 0, tol_val)
names(ret) = as.character(lam2s)
ret = rbindlist(ret, use.names = T, idcol = "lam2")
ret[,group:="all"]
ggplot(ret)+geom_raster(aes(bin1,bin2,fill=pmin(phihat,3)))+geom_raster(aes(bin2,bin1,fill=pmin(signal,3)))+facet_wrap(~factor(lam2))+coord_fixed()+scale_fill_gradient2()
ggplot(ret)+geom_raster(aes(bin1,bin2,fill=pmin(signal,3)),data=ret.old)+geom_raster(aes(bin2,bin1,fill=pmin(signal,3)))+facet_wrap(~factor(lam2))+coord_fixed()+scale_fill_gradient2()
#
retg1=binless:::fast_binless_eval_cv(out, lam2s, 1, tol_val)
names(retg1) = as.character(lam2s)
retg1 = rbindlist(retg1, use.names = T, idcol = "lam2")
retg1[,group:="g1"]
#
retg2=binless:::fast_binless_eval_cv(out, lam2s, 2, tol_val)
names(retg2) = as.character(lam2s)
retg2 = rbindlist(retg2, use.names = T, idcol = "lam2")
retg2[,group:="g2"]
#
retg=rbind(retg1,retg2)[phihat==-100]
setkey(retg,lam2,name,bin1,bin2)
retall=dcast(rbindlist(list(all=ret,cv=retg),use.names = T, idcol = "mode"), lam2+name+bin1+bin2+observed+biases+decay~mode,
             value.var=c("expected","signal","binless","group","phihat","weights"))
#retall[,lam2:=ordered(lam2,lam2s)]
retall[,lam2:=as.numeric(lam2)]
retall[,sqerr_cv:=(phihat_all-signal_cv)**2]
retall[,sqerr_all:=(phihat_all-signal_all)**2]
scores = retall[,.(cv_group=mean(weights_all*sqerr_cv),n_group=.N),by=c("lam2","group_cv")][
  ,.(cv=sum(n_group*cv_group)/sum(n_group), cv.sd=sum(n_group*cv_group**2)/sum(n_group)),keyby=lam2]
scores[,cv.sd:=sqrt(cv.sd-cv**2)]
ggplot(scores)+geom_line(aes(lam2,cv))+geom_errorbar(aes(lam2,ymin=cv-cv.sd,ymax=cv+cv.sd))
#
ggplot(retall[,.(CV=sum(weights_all*sqerr_cv),sumsq=sum(weights_all*sqerr_all)),by=lam2])+
  geom_line(aes(lam2,CV,colour="CV",group=1))+geom_line(aes(lam2,sumsq,colour="sumsq",group=1))
#
ggplot(retall[,.(CV=sum(pmin(weights_all,10)*sqerr_cv),sumsq=sum(pmin(weights_all,10)*sqerr_all)),by=lam2])+
  geom_line(aes(lam2,CV,colour="CV",group=1))+geom_line(aes(lam2,sumsq,colour="sumsq",group=1))
#
ggplot(retall[,.(CV=sum(weights_all*sqerr_cv/exp(decay)),sumsq=sum(weights_all*sqerr_all/exp(decay))),by=lam2])+
  geom_line(aes(lam2,CV,colour="CV",group=1))+geom_line(aes(lam2,sumsq,colour="sumsq",group=1))
#
ggplot(retall[,.(CV=sum(sqerr_cv),sumsq=sum(sqerr_all)),by=lam2])+
  geom_line(aes(lam2,CV,colour="CV",group=1))+geom_line(aes(lam2,sumsq,colour="sumsq",group=1))
#
ggplot(retall[,.(msq_wt=pmin(mean(weights_all*sqerr_cv),100),msq_flat=pmin(mean(sqerr_cv),100)),by=c("bin1","bin2","good")])+geom_raster(aes(bin1,bin2,fill=msq_wt))+
  geom_raster(aes(bin2,bin1,fill=msq_flat))+coord_fixed()+scale_fill_gradient2()+facet_wrap(~good)
#
ggplot(ret)+geom_raster(aes(bin1,bin2,fill=pmin(phihat,max(signal))))+geom_raster(aes(bin2,bin1,fill=signal))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()
ggplot(ret)+geom_raster(aes(bin1,bin2,fill=observed))+geom_raster(aes(bin2,bin1,fill=binless))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()
ggplot(ret[lam2==lam2[.N]])+geom_raster(aes(bin1,bin2,fill=expected))+geom_raster(aes(bin2,bin1,fill=weights))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()

ggplot(retall[bin2<20])+geom_raster(aes(bin1,bin2,fill=group_cv))+geom_raster(aes(bin2,bin1,fill=group_all))+facet_wrap(~lam2)+coord_fixed()
ggplot(retall[bin2<20])+geom_raster(aes(bin1,bin2,fill=signal_cv))+geom_raster(aes(bin2,bin1,fill=signal_all))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()
ggplot(retall)+geom_raster(aes(bin1,bin2,fill=signal_cv))+geom_raster(aes(bin2,bin1,fill=signal_all))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()
ggplot(retall[bin1>225&bin1<250&bin2>425&bin2<450])+geom_raster(aes(bin1,bin2-150,fill=signal_cv))+geom_raster(aes(bin2-150,bin1,fill=signal_all))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()
ggplot(retall[bin1>225&bin1<250&bin2>425&bin2<450])+geom_raster(aes(bin1,bin2-150,fill=phihat))+geom_raster(aes(bin2-150,bin1,fill=signal_all))+facet_wrap(~lam2)+coord_fixed()+scale_fill_gradient2()

