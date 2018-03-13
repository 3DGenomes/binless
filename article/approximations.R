library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)


### evaluate the impact of approximating the zeros means in optimized binless iterations

# adapted from gauss_common_muhat_mean in common.R
predict_all = function(cs, counts) {
  sbins=cs@settings$sbins
  init=cs@par
  ### positive counts (twice, in both directions)
  #compute means
  cpos=copy(counts)
  bsub=cs@biases[,.(id)]
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  cpos=merge(bsub[,.(id1=id,log_iota,log_rho)],cpos,by="id1",all.x=F,all.y=T)
  cpos=merge(bsub[,.(id2=id,log_iota,log_rho)],cpos,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  cpos=merge(cbind(cs@design[,.(name)],eC=init$eC), cpos, by="name",all.x=F,all.y=T)
  cpos[,c("bin1","bin2","dbin"):=
         list(cut(pos1, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(pos2, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cpos=merge(cpos,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  cpos[,log_mu.base:=eC + log_decay]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  cpos[,log_mu.base:=NULL]
  cpos=rbind(cpos[,.(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd", id2=id2, variable="close",
                                    count=contact.close, lmu.nosig=lmu.close, nobs=1, eC, log_decay, log_bias=log_rho1)],
             cpos[,.(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd", id2=id2, variable="far",
                                    count=contact.far,   lmu.nosig=lmu.far,   nobs=1, eC, log_decay, log_bias=log_iota1)],
             cpos[,.(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd", id2=id2, variable="down",
                                    count=contact.down,  lmu.nosig=lmu.down,  nobs=1, eC, log_decay, log_bias=log_rho1)],
             cpos[,.(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd", id2=id2, variable="up",
                                    count=contact.up,    lmu.nosig=lmu.up,    nobs=1, eC, log_decay, log_bias=log_iota1)],
             cpos[,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev", id2=id1, variable="far",
                                    count=contact.far,   lmu.nosig=lmu.far,   nobs=1, eC, log_decay, log_bias=log_rho2)],
             cpos[,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev", id2=id1, variable="close",
                                    count=contact.close, lmu.nosig=lmu.close, nobs=1, eC, log_decay, log_bias=log_iota2)],
             cpos[,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev", id2=id1, variable="down",
                                    count=contact.down,  lmu.nosig=lmu.down,  nobs=1, eC, log_decay, log_bias=log_rho2)],
             cpos[,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev", id2=id1, variable="up",
                                    count=contact.up,    lmu.nosig=lmu.up,    nobs=1, eC, log_decay, log_bias=log_iota2)])
  cts=cpos
  setkeyv(cts,key(cs@zeros))
  ### add signal
  signal = binless:::get_signal_matrix(cs, resolution = sbins[2]-sbins[1], groups=cs@experiments[,.(name,groupname=name)])
  signal=rbind(signal[,.(name,bin1,bin2,phi)],signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
  cts=signal[cts,,on=c("name","bin1","bin2")]
  cts[,mu:=exp(lmu.nosig+phi)]
  ### finalize
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  setkey(cts,name,id1,pos1,bin1,bin2,dbin,dir,cat)
  #ggplot(cts[name=="T47D es 60 MboI 1"&cat=="contact R"])+geom_line(aes(dbin,log_decay,colour=count>0,group=count>0))
  return(cts)
}

#load dataset and create CS object
load("data/rao_HiCall_GM12878_SELP_150k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")

#generate dense counts matrix with all observations, including zeros
counts=binless:::subsample_counts(cs, cs@biases[,.N*(.N-1)/2])

#generate some statistics on the dataset
counts_stats=melt(counts,id.vars = c("name","id1","id2","pos1","pos2","distance"))[,.(ncounts=.N),keyby=c("variable","value")]
counts_stats=counts_stats[CJ(variable=variable,value=0:max(value),unique=T)]
counts_stats[,value:=factor(value)]
counts_stats[is.na(ncounts),ncounts:=0]

counts_stats[value==0,sum(ncounts)]/counts_stats[,sum(ncounts)]*100
counts_stats[value==1,sum(ncounts)]/counts_stats[,sum(ncounts)]*100
ggplot(counts_stats)+geom_bar(aes(value,ncounts,fill=variable),position="dodge",stat="identity")+scale_y_log10()

#normalize using optimized binless
cs <- normalize_binless(cs, ngibbs = 15, ncores = 1, base.res = 5000, bg.steps = 5, tol = 1e-1)

#get exact and approximate predictions, merge to each count
pred.exact=predict_all(cs,counts)
pred.approx=binless:::gauss_common_muhat_mean(cs,cs@zeros,cs@settings$sbins)
setkeyv(pred.approx,c(key(pred.approx),"count"))
setkeyv(pred.exact,key(pred.approx))
pred=merge(pred.approx,pred.exact,all=T,suffixes=c(".approx",".exact"))

#compare means
pred[,summary(eC.exact-eC.approx)]
pred[,summary(log_decay.exact-log_decay.approx)]
pred[,summary(log_bias.exact-log_bias.approx)]
pred[,summary(phi.exact-phi.approx)]

pred[,summary(lmu.nosig.exact-lmu.nosig.approx)]
pred[,summary(lm(lmu.nosig.approx~lmu.nosig.exact+0))]
ggplot(pred)+geom_point(aes(lmu.nosig.exact,lmu.nosig.approx),alpha=0.01)+stat_function(fun=identity)
ggplot(aes(lmu.nosig.exact,lmu.nosig.approx),data=pred)+geom_density_2d()+ stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+stat_function(fun=identity)

#show approximation improves by averaging
pred.avg=pred[,.(lmu.nosig.exact=mean(lmu.nosig.exact),lmu.nosig.approx=mean(lmu.nosig.approx)),by=c("name","id1","bin2","dbin")]
pred.avg[,summary(lmu.nosig.exact-lmu.nosig.approx)]
pred.avg[,summary(lm(lmu.nosig.approx~lmu.nosig.exact+0))]
ggplot(pred.avg)+geom_point(aes(lmu.nosig.exact,lmu.nosig.approx),alpha=0.1)+stat_function(fun=identity)
ggplot(aes(lmu.nosig.exact,lmu.nosig.approx),data=pred.avg)+geom_density_2d()+ stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+stat_function(fun=identity)

pred.avg2=pred[,.(lmu.nosig.exact=mean(lmu.nosig.exact),lmu.nosig.approx=mean(lmu.nosig.approx)),by=c("name","id1","bin2")]
pred.avg2[,summary(lmu.nosig.exact-lmu.nosig.approx)]
pred.avg2[,summary(lm(lmu.nosig.approx~lmu.nosig.exact+0))]
ggplot(pred.avg2)+geom_point(aes(lmu.nosig.exact,lmu.nosig.approx),alpha=0.1)+stat_function(fun=identity)
ggplot(aes(lmu.nosig.exact,lmu.nosig.approx),data=pred.avg2)+geom_density_2d()+ stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+stat_function(fun=identity)


### compare optimized and slow binless

#First iteration
base.res=5000
tol_val=1e-1
bg_steps=5
cs2 <- normalize_binless(cs, ngibbs = 0, ncores = 1, base.res = base.res, bg.steps = bg_steps, tol = tol_val)
cts.common = binless:::gauss_common_muhat_mean(cs2, cs2@zeros, cs2@settings$sbins)
cts.common = cts.common[,.(name,bin1=pmin(bin1,bin2),bin2=pmax(bin1,bin2),count,z,var,mu,lmu.nosig,
                    log_decay,log_bias,nobs=nobs/2)]
pred.opt=cts.common[,.(observed=sum(count*nobs),nobs=sum(nobs),log_background=weighted.mean(lmu.nosig,nobs),weight=sum(nobs/var)),by=c("name","bin1","bin2")]
pred.opt[,c("name","bin1","bin2"):=list(unclass(name),unclass(bin1),unclass(bin2))]
pred.opt=pred.opt[bin1<max(bin2)&bin2<max(bin2)] #in fast binless, last row is discarded
mat = binless:::bin_data(cs,resolution=base.res)
nouter=0
lam2=5
out=binless:::fast_binless(mat, mat[,nlevels(bin1)], lam2, cs2@par$alpha, nouter, tol_val, bg_steps)
pred.fast=as.data.table(out$mat)[,.(name,bin1,bin2,distance,observed,nobs,log_background)]
pred.fast[,weight:=nobs/(1/exp(log_background)+1/cs2@par$alpha)]

#compare observed
pred=merge(pred.opt,pred.fast,all=T,by=c("name","bin1","bin2"),suffixes=c(".opt",".fast"))
pred[,summary(observed.opt-observed.fast)]
ggplot(pred)+geom_raster(aes(bin1,bin2,fill=observed.opt-observed.fast))+scale_fill_gradient2() #all equal except most distant corner, on purpose
#compare nobs
pred[,summary(nobs.opt-nobs.fast)]
ggplot(pred)+geom_raster(aes(bin1,bin2,fill=nobs.opt-nobs.fast))+scale_fill_gradient2() #all equal except most distant corner, on purpose
ggplot(pred)+geom_point(aes(nobs.opt,nobs.fast,colour=distance<=base.res/2+cs2@settings$dmin))+stat_function(fun=identity)
#compare log_background
pred[,summary(log_background.opt-log_background.fast)]
ggplot(pred)+geom_raster(aes(bin1,bin2,fill=log_background.opt-log_background.fast))+scale_fill_gradient2() #all equal except most distant corner, on purpose
ggplot(pred)+geom_point(aes(log_background.opt,log_background.fast,colour=distance<=base.res/2+cs2@settings$dmin))+stat_function(fun=identity)
pred[,summary(lm(log_background.fast~log_background.opt+0))]
#compare weights
pred[,summary((weight.opt-weight.fast)/weight.opt)]
ggplot(pred)+geom_raster(aes(bin1,bin2,fill=(weight.opt-weight.fast)/weight.opt))+scale_fill_gradient2() #all equal except most distant corner, on purpose
ggplot(pred)+geom_point(aes(weight.opt,weight.fast,colour=distance<=base.res/2+cs2@settings$dmin))+stat_function(fun=identity)+scale_x_log10()+scale_y_log10()
pred[,summary(lm(log_background.fast~log_background.opt+0))]




