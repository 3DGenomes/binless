library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)
library(scam)
library(Hmisc)
library(foreach)
library(doParallel)
registerDoParallel()

setwd("/home/yannick/simulations/spline_stan")

stan_matrix_to_datatable = function(opt, x) {
  vals=data.table(opt)
  vals[,x:=x]
  melt(data.table(vals), id.vars="x")
}

optimize_all_meanfield = function(model, biases, counts, meanfield, maxcount, Krow=1000, Kdiag=10, 
                               lambda_nu=1, lambda_delta=1, lambda_diag=1, iter=10000, verbose=T) {
  counts_ex=counts[count>maxcount]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  optimizing(model, data = list( Krow=Krow, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                 danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                 Kdiag=Kdiag, Nexpl=counts[,.N],
                                 counts=t(data.matrix(counts_ex[,.(contact.close,contact.far,contact.up,contact.down)])),
                                 cidx=t(data.matrix(counts_ex[,.(id1,id2)])),
                                 Nl=mf$Nkl[,.N], Nkl_count=mf$Nkl[,count], Nkl_cidx=mf$Nkl[,id], Nkl_N=mf$Nkl[,N],
                                 Nr=mf$Nkr[,.N], Nkr_count=mf$Nkr[,count], Nkr_cidx=mf$Nkr[,id], Nkr_N=mf$Nkr[,N],
                                 Nd=mf$Nkd[,.N], Nkd_count=mf$Nkd[,count], Nkd_d=mf$Nkl[,mdist], Nkd_N=mf$Nkd[,N],
                                 lambda_nu=lambda_nu, lambda_delta=lambda_delta, lambda_diag=lambda_diag),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

predict_full = function(model, biases, counts, opt, Kdiag=10, verbose=T) {
  optimizing(model, data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], N=counts[,.N],
                                 counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                 cidx=t(data.matrix(counts[,.(id1,id2)])),
                                 eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
                                 beta_diag=opt$par$beta_diag, alpha=opt$par$alpha),
             as_vector=F, hessian=F, iter=1, verbose=verbose)
}


bin_counts = function(counts, biases, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL, normalized=T) {
  if (is.null(b1)) b1=counts[,min(pos1)]-1
  if (is.null(b2)) b2=counts[,min(pos2)]-1
  if (is.null(e1)) e1=counts[,max(pos1)]+1
  if (is.null(e2)) e2=counts[,max(pos2)]+1
  bins1=seq(b1,e1+resolution,resolution)
  bins2=seq(b2,e2+resolution,resolution)
  #
  counts.wrapped=rbindlist(list(counts[,.(pos1,pos2,contact.close,log_decay,log_mean_cclose)],
                                counts[,.(pos1,pos2,contact.far,log_decay,log_mean_cfar)],
                                counts[,.(pos1,pos2,contact.up,log_decay,log_mean_cup)],
                                counts[,.(pos1,pos2,contact.down,log_decay,log_mean_cdown)]))
  setnames(counts.wrapped, c("pos1","pos2","count","log_decay","log_mean"))
  if (normalized==T) {
    sub = counts.wrapped[,.(pos1,pos2,bin1=cut2(pos1, bins1, oneval=F, onlycuts=T, digits=10, minmax=F),
                            bin2=cut2(pos2,bins2, oneval=F, onlycuts=T, digits=10, minmax=F),
                            weight=exp(log_mean-log_decay))
                         ][,.(N=sum(weight, na.rm = T)),by=c("bin1","bin2")]
  } else {
    sub = counts.wrapped[,.(pos1,pos2,bin1=cut2(pos1, bins1, oneval=F, onlycuts=T, digits=10, minmax=F),
                            bin2=cut2(pos2,bins2, oneval=F, onlycuts=T, digits=10, minmax=F),
                            weight=count)
                         ][,.(N=sum(weight, na.rm = T)),by=c("bin1","bin2")]
  }
  #remove NAs in case a zoom was performed
  sub=sub[complete.cases(sub)]
  #divide by number of rsites in each bin
  if (normalized==T) {
    ns1 = biases[,.(bin1=cut2(pos, bins1, oneval=F, onlycuts=F, digits=10))][,.(nsites1=.N),keyby=bin1]
    setkey(sub, bin1)
    sub=ns1[sub]
    ns2 = biases[,.(bin2=cut2(pos, bins2, oneval=F, onlycuts=F, digits=10))][,.(nsites2=.N),keyby=bin2]
    setkey(sub, bin2)
    sub=ns2[sub]
    sub[,N:=N/(nsites1*nsites2)]
  }
  #write begins/ends
  bin1.begin=sub[,bin1]
  bin1.end=sub[,bin1]
  bin2.begin=sub[,bin2]
  bin2.end=sub[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[[,]")[2][[1]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[[,)]")[2][[1]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[[,]")[2][[1]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[[,)]")[2][[1]]
  sub[,begin1:=as.integer(as.character(bin1.begin))]
  sub[,end1:=as.integer(as.character(bin1.end))]
  sub[,begin2:=as.integer(as.character(bin2.begin))]
  sub[,end2:=as.integer(as.character(bin2.end))]
  return(sub)
}

bin_for_mean_field = function(biases, counts, distance_bins_per_decade=10) {
  stopifnot(counts[id1>=id2,.N]==0)
  mcounts=melt(counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0]
  ### accumulate counts
  ci=mcounts[,.(id=id1,count,category)][,.N,by=c("id","count","category")]
  cj=mcounts[,.(id=id2,count,category)][,.N,by=c("id","count","category")]
  nsites=biases[,.N]
  ### make histograms for biases
  Nkl=dcast(rbind(ci[category=="contact.up",.(id,count,category="Ni.up",N)],
                  ci[category=="contact.far",.(id,count,category="Ni.far",N)],
                  cj[category=="contact.up",.(id,count,category="Nj.up",N)],
                  cj[category=="contact.close",.(id,count,category="Nj.close",N)]),
            ...~category, value.var="N", fill=0)[,.(id,count,N=Ni.far+Ni.up+Nj.up+Nj.close)]
  Nkl=rbind(Nkl,Nkl[,.(count=0,N=2*nsites-sum(N)),by=id]) #each rsite is counted twice
  setkey(Nkl,N,id,count)
  Nkr=dcast(rbind(ci[category=="contact.close",.(id,count,category="Ni.close",N)],
                  ci[category=="contact.down",.(id,count,category="Ni.down",N)],
                  cj[category=="contact.far",.(id,count,category="Nj.far",N)],
                  cj[category=="contact.down",.(id,count,category="Nj.down",N)]),
            ...~category, value.var="N", fill=0)[,.(id,count,N=Ni.close+Ni.down+Nj.far+Nj.down)]
  Nkr=rbind(Nkr,Nkr[,.(count=0,N=2*nsites-sum(N)),by=id]) #each rsite is counted twice
  setkey(Nkr,N,id,count)
  ### make histogram for distance
  #make distance bins and their factor
  stepsz=1/distance_bins_per_decade
  dbins=10**seq(0,biases[,log10(max(pos)-min(pos))]+stepsz,stepsz)
  mcounts[,distance:=abs(pos1-pos2)]
  mcounts[,bdist:=cut(distance,dbins,ordered_result=T,right=F)]
  #Count positive counts in these bins
  Nkd = mcounts[,.N,keyby=c("bdist","count")]
  Nkd[,mdist:=sqrt(dbins[unclass(bdist)+1]*dbins[unclass(bdist)])]
  #Count the number of crossings per distance bin
  positions=biases[,pos]
  npos=length(positions)
  ncrossings <- rowSums(sapply(1:(npos-1), #this loop is still 5x faster in python
                               function(i){hist(positions[(i+1):npos]-positions[i],breaks=dbins,plot=F)$counts}
  ))
  #deduce zero counts
  Nkz = data.table(bdist=Nkd[,ordered(levels(bdist), levels(bdist))], ncrossings=ncrossings, key="bdist")
  Nkz = Nkd[,.(nnz=sum(N)),by=bdist][Nkz[ncrossings>0]]
  Nkz[,mdist:=sqrt(dbins[unclass(bdist)+1]*dbins[unclass(bdist)])]
  Nkz[is.na(nnz),nnz:=0]
  Nkd = rbind(Nkd, Nkz[,.(bdist,mdist,count=0,N=4*ncrossings-nnz)]) # one crossing for each of 4 count types
  setkey(Nkd,N,bdist,count)
  stopifnot(Nkd[count==0,all(N>=0)])
  #plot decay
  #ggplot(Nkd[,.(decay=sum(N*count)/sum(N)),by=bdist])+geom_point(aes(bdist,decay))+scale_y_log10()+
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(list(Nkd=Nkd, Nkr=Nkr, Nkl=Nkl))
}

gammaQuery=function(mu,theta,count){
  alpha1=theta
  alpha2=theta+count
  beta1=theta/mu
  beta2=beta1+1
  return(function(x){dgamma(x,alpha2,rate=beta2)*pgamma(x,alpha1,rate=beta1)})
}
igammaQuery_c=function(mu,theta,threshold=0.95){
  function(x){integrate(gammaQuery(mu=mu,theta=theta,count=x),0,+Inf)$value-threshold}
}
#uniroot(igammaQuery_c(mu=2,theta=2,threshold=0.9),c(0,10))$root
#uniroot(igammaQuery_c(mu=2,theta=20,threshold=0.9),c(10,30))$root
igammaQuery_theta=function(mu,count,threshold=0.95){
  function(x){integrate(gammaQuery(mu=mu,theta=x,count=count),0,+Inf)$value-threshold}
}
#uniroot(igammaQuery_theta(mu=2,count=3,threshold=0.9),c(0,30))$root
#uniroot(igammaQuery_theta(mu=2,count=10,threshold=0.9),c(1,30))$root



biases=fread("data/rao_HICall_chr20_all_biases.dat") #157438
setkey(biases,id)
counts=fread("data/rao_HICall_chr20_all_counts.dat")
both=fread("data/rao_HICall_chr20_all_both.dat")


biases=fread("data/rao_HIC035_chr20_all_biases.dat") #16345
setkey(biases,id)
counts=fread("data/rao_HIC035_chr20_all_counts.dat")
load("data/rao_HIC035_chr20_all_meanfield.RData")
meanfield = bin_for_mean_field(biases, counts, distance_bins_per_decade = 10)
save(meanfield, file = "data/rao_HIC035_chr20_all_meanfield.RData")
both=fread("data/rao_HIC035_chr20_all_both.dat")


biases=fread("data/rao_HICall_chrX_73780165-74230165_biases.dat") #1201
setkey(biases,id)
counts=fread("data/rao_HICall_chrX_73780165-74230165_counts.dat")
load("data/rao_HICall_chrX_73780165-74230165_meanfield.RData")
#meanfield = bin_for_mean_field(biases, counts, distance_bins_per_decade = 10)
#save(meanfield, file = "data/rao_HICall_chrX_73780165-74230165_meanfield.RData")
both=fread("data/rao_HICall_chrX_73780165-74230165_both.dat")

biases=fread("data/rao_HICall_chr20_35000000-36000000_biases.dat") #3215
setkey(biases,id)
counts=fread("data/rao_HICall_chr20_35000000-36000000_counts.dat")
both=fread("data/rao_HICall_chr20_35000000-36000000_both.dat")

biases=fread("data/rao_HIC035_chr20_35000000-36000000_biases.dat") #212
setkey(biases,id)
counts=fread("data/rao_HIC035_chr20_35000000-36000000_counts.dat")
both=fread("data/rao_HIC035_chr20_35000000-36000000_both.dat")


biases=fread("data/caulo_all_biases.dat") #2010
setkey(biases,id)
counts=fread("data/caulo_all_counts.dat")
both=fread("data/caulo_all_both.dat")
counts=counts[sample(.N,min(25000,.N))]

biases=fread("data/caulo_3000000-4000000_biases.dat") #525
setkey(biases,id)
counts=fread("data/caulo_3000000-4000000_counts.dat")
both=fread("data/caulo_3000000-4000000_both.dat")




#### optimization wihout prior guesses
counts.sub=counts[sample(.N,min(.N,100000))]
smfit = stan_model(file = "sparse_cs_norm_fit_meanfield.stan")
system.time(op <- optimize_all_meanfield(smfit, biases, counts.sub, Krow=4000, Kdiag=10,
                                      lambda_nu=10, lambda_delta=10, lambda_diag=.1, verbose = T))
save(op, file = "data/rao_HICall_chrX_73780165-74230165_op_lambda10.RData")
#compare deviances
c(100*(fit$null.deviance-fit$deviance)/fit$null.deviance, op$par$deviance_proportion_explained)
#compare dispersions
c(fit$family$getTheta(T), op$par$alpha)
#compare decay
a=cbind(counts.sub,data.table(decay=exp(op$par$log_decay),
                              cclose=exp(op$par$log_mean_cclose),
                              cfar=exp(op$par$log_mean_cfar),
                              cup=exp(op$par$log_mean_cup),
                              cdown=exp(op$par$log_mean_cdown)))
ggplot(a[pos2-pos1>1e4][sample(.N,min(.N,10000))])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(pos2-pos1, contact.close*decay/cclose), colour="blue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.far*decay/cclose), colour="green", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.up*decay/cclose), colour="darkblue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.down*decay/cclose), colour="darkgreen", alpha=0.01)+
  geom_line(aes(pos2-pos1, decay) ,colour="pink")
#nu
a=cbind(biases,data.table(opt=exp(op$par$log_nu)))
#pbegin=35100000 
#pend=35200000 
#pbegin=3166716
#pend=3191637
pbegin=73780165
pend=74230165
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op$par$eDE)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op$par$eDE)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op$par$eRJ)),colour="red")+
  geom_point(aes(pos, opt),colour="blue")
#delta
a=cbind(biases,data.table(opt=exp(op$par$log_delta),
                          log_nu=op$par$log_nu))
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/(exp(log_nu+op$par$eDE))),colour="orange")+
  geom_point(aes(pos, dangling.R/(exp(log_nu+op$par$eDE))),colour="pink")+
  geom_point(aes(pos, rejoined/(exp(log_nu+op$par$eRJ))),colour="red")+
  geom_point(aes(pos, opt),colour="orange", shape=0)+
  geom_point(aes(pos, 1/opt),colour="pink", shape=0)
#all 3 types
a=cbind(biases,data.table(RJ=exp(op$par$eRJ+op$par$log_nu),
                          DL=exp(op$par$eDE+op$par$log_nu+op$par$log_delta),
                          DR=exp(op$par$eDE+op$par$log_nu-op$par$log_delta)))
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, rejoined), colour="red")+ geom_point(aes(pos, RJ), shape=0, colour="red")+
  geom_point(aes(pos, dangling.L), colour="orange")+ geom_point(aes(pos, DL), shape=0, colour="orange")+
  geom_point(aes(pos, dangling.R), colour="pink")+ geom_point(aes(pos, DR), shape=0, colour="pink")

## predict on full dataset
smpred = stan_model(file = "sparse_cs_norm_predict.stan")
system.time(op.pred <- predict_full(smpred, biases, counts, op, Kdiag=10, verbose = T))
counts[,log_decay:=op.pred$par$log_decay]
counts[,log_mean_cup:=op.pred$par$log_mean_cup]
counts[,log_mean_cdown:=op.pred$par$log_mean_cdown]
counts[,log_mean_cfar:=op.pred$par$log_mean_cfar]
counts[,log_mean_cclose:=op.pred$par$log_mean_cclose]
write.table(counts, file = "data/rao_HICall_chrX_73780165-74230165_counts_fitted_disp3.84_dev_14.3_lambda10.dat", quote = F, row.names = F)
#divide by rsites only
counts[,log_decay:=0]
counts[,log_mean_cup:=log(contact.up)]
counts[,log_mean_cdown:=log(contact.down)]
counts[,log_mean_cfar:=log(contact.far)]
counts[,log_mean_cclose:=log(contact.close)]


#visualize binned matrices
binned = bin_counts(counts=counts, biases=biases, resolution=5000, normalized=F)
ggplot(binned, aes(begin1,begin2, fill=log(N)))+geom_raster()+scale_fill_gradient(low="white", high="black")
ggsave(filename = "images/rao_HICall_chrX_73780165-74230165_5k_raw.png", width=10, height=7.5)
write.table(binned[,.(begin1,end1,begin2,end2,N)], file = "data/binned_Rao_MboI_chrX_73780165-74230165_raw_5k.dat", quote = F, row.names = F)
#
binned = bin_counts(counts=counts, biases=biases, resolution=5000, normalized=T)
ggplot(binned, aes(begin1,begin2, fill=log(N)))+geom_raster()+scale_fill_gradient(low="white", high="black")
ggsave(filename = "images/rao_HICall_chrX_73780165-74230165_5k_div_nsites.png", width=10, height=7.5)
write.table(binned[,.(begin1,end1,begin2,end2,N)], file = "data/binned_Rao_MboI_chrX_73780165-74230165_div_nsites.dat", quote = F, row.names = F)



