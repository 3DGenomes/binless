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
                               lambda_nu=1, lambda_delta=1, lambda_diag=1, iter=10000, verbose=T, mincount=-1) {
  cclose=counts[contact.close>maxcount,.(id1,id2,count=contact.close)][count>mincount]
  cfar=counts[contact.far>maxcount,.(id1,id2,count=contact.far)][count>mincount]
  cup=counts[contact.up>maxcount,.(id1,id2,count=contact.up)][count>mincount]
  cdown=counts[contact.down>maxcount,.(id1,id2,count=contact.down)][count>mincount]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])),
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])),
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])),
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])),
               Nl=mf$Nkl[,.N], Nkl_count=mf$Nkl[,count], Nkl_cidx=mf$Nkl[,id], Nkl_N=mf$Nkl[,N], Nkl_levels=mf$Nkl[,sum(diff(N)!=0)+1],
               Nr=mf$Nkr[,.N], Nkr_count=mf$Nkr[,count], Nkr_cidx=mf$Nkr[,id], Nkr_N=mf$Nkr[,N], Nkr_levels=mf$Nkr[,sum(diff(N)!=0)+1],
               Nd=mf$Nkd[,.N], Nkd_count=mf$Nkd[,count], Nkd_d=mf$Nkd[,mdist], Nkd_N=mf$Nkd[,N], Nkd_levels=mf$Nkd[,sum(diff(N)!=0)+1],
               lambda_nu=lambda_nu, lambda_delta=lambda_delta, lambda_diag=lambda_diag)
  if (data$Nl==0) data$Nkl_levels=0
  if (data$Nr==0) data$Nkr_levels=0
  if (data$Nd==0) data$Nkd_levels=0
  message("Biases      : ", biases[,.N])
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  message("Left counts : ", mf$Nkl[,.N], " (", data$Nkl_levels, " levels)")
  message("Right counts: ", mf$Nkr[,.N], " (", data$Nkr_levels, " levels)")
  message("Decay counts: ", mf$Nkd[,.N], " (", data$Nkd_levels, " levels)")
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

predict_full = function(model, biases, counts, opt, Kdiag=10, verbose=T) {
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], N=counts[,.N],
               counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
               cidx=t(data.matrix(counts[,.(id1,id2)])),
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag=opt$par$beta_diag, alpha=opt$par$alpha)
  optimizing(model, data = data, as_vector=F, hessian=F, iter=1, verbose=verbose)
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

dset_statistics = function(biases,counts){
  message("Mean counts")
  message("   Rejoined  : ", biases[,mean(rejoined)])
  message("   Dangling L: ", biases[,mean(dangling.L)])
  message("   Dangling R: ", biases[,mean(dangling.R)])
  message("   C. close  : ", counts[,mean(contact.close)])
  message("   C. far    : ", counts[,mean(contact.far)])
  message("   C. up     : ", counts[,mean(contact.up)])
  message("   C. down   : ", counts[,mean(contact.down)])
  message("Median counts")
  message("   Rejoined  : ", biases[,median(rejoined)])
  message("   Dangling L: ", biases[,median(dangling.L)])
  message("   Dangling R: ", biases[,median(dangling.R)])
  message("   C. close  : ", counts[,median(contact.close)])
  message("   C. far    : ", counts[,median(contact.far)])
  message("   C. up     : ", counts[,median(contact.up)])
  message("   C. down   : ", counts[,median(contact.down)])
  message("Percent of zero counts")
  message("   Rejoined  : ", biases[rejoined==0,.N]/biases[,.N]*100)
  message("   Dangling L: ", biases[dangling.L==0,.N]/biases[,.N]*100)
  message("   Dangling R: ", biases[dangling.R==0,.N]/biases[,.N]*100)
  message("   C. close  : ", counts[contact.close==0,.N]/counts[,.N]*100)
  message("   C. far    : ", counts[contact.far==0,.N]/counts[,.N]*100)
  message("   C. up     : ", counts[contact.up==0,.N]/counts[,.N]*100)
  message("   C. down   : ", counts[contact.down==0,.N]/counts[,.N]*100)
}

generate_fake_dataset = function(num_rsites=10, genome_size=10000, eC=.1, eRJ=.4, eDE=.8, alpha=10) {
  #place rsites
  biases=data.table(id=seq(num_rsites),
                    pos=cumsum(rmultinom(n=1, size=genome_size, prob=rep(1,num_rsites+1)))[1:num_rsites])
  setkey(biases,id)
  #build biases
  biases[,true_log_nu:=sin(pos/1000)+(pos-genome_size/2)/genome_size]
  biases[,true_log_nu:=true_log_nu-mean(true_log_nu)]
  ggplot(biases,aes(pos,true_log_nu))+geom_point()+geom_line()
  biases[,true_log_delta:=-(pos-genome_size/2)/genome_size+sin(pos/2100)]
  biases[,true_log_delta:=true_log_delta-mean(true_log_delta)]
  ggplot(biases,aes(pos,true_log_delta))+geom_point()+geom_line()
  biases[,true_log_mean_RJ:=eRJ+true_log_nu]
  biases[,true_log_mean_DL:=eDE+true_log_nu+true_log_delta]
  biases[,true_log_mean_DR:=eDE+true_log_nu-true_log_delta]
  #draw dangling/rejoined
  biases[,dangling.L:=rnbinom(.N, mu=exp(true_log_mean_DL), size=alpha)]
  biases[,dangling.R:=rnbinom(.N, mu=exp(true_log_mean_DR), size=alpha)]
  biases[,rejoined:=rnbinom(.N, mu=exp(true_log_mean_RJ), size=alpha)]
  #report rsites in counts
  counts=CJ(biases[,paste(id,pos,true_log_nu,true_log_delta)],biases[,paste(id,pos,true_log_nu,true_log_delta)])
  counts[,c("id1","pos1","true_log_nu1","true_log_delta1"):=tstrsplit(V1, " ")]
  counts[,c("id2","pos2","true_log_nu2","true_log_delta2"):=tstrsplit(V2, " ")]
  counts[,c("id1","id2","pos1","pos2","V1","V2"):=list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
  counts[,c("true_log_nu1","true_log_nu2","true_log_delta1","true_log_delta2"):=list(
    as.numeric(true_log_nu1), as.numeric(true_log_nu2), as.numeric(true_log_delta1), as.numeric(true_log_delta2))]
  counts=counts[pos1<pos2]
  setkey(counts, id1, id2)
  #build decay
  counts[,distance:=pos2-pos1]
  counts[,true_log_decay:=5*dnorm(log(distance), mean=log(min(distance)), sd=log(max(distance))/10)-0.2*log(distance)]
  counts[,true_log_decay:=true_log_decay-mean(true_log_decay)]
  ggplot(dset$counts)+geom_point(aes(distance,exp(true_log_decay)))+scale_x_log10()+scale_y_log10()
  counts[,base_count:=eC+true_log_decay+true_log_nu1+true_log_nu2]
  counts[,true_log_mean_cclose:=base_count-true_log_delta1+true_log_delta2]
  counts[,true_log_mean_cfar:=base_count+true_log_delta1-true_log_delta2]
  counts[,true_log_mean_cup:=base_count+true_log_delta1+true_log_delta2]
  counts[,true_log_mean_cdown:=base_count-true_log_delta1-true_log_delta2]
  #draw counts
  counts[,contact.close:=rnbinom(.N, mu=exp(true_log_mean_cclose), size=alpha)]
  counts[,contact.far:=rnbinom(.N, mu=exp(true_log_mean_cfar), size=alpha)]
  counts[,contact.up:=rnbinom(.N, mu=exp(true_log_mean_cup), size=alpha)]
  counts[,contact.down:=rnbinom(.N, mu=exp(true_log_mean_cdown), size=alpha)]
  ggplot(dset$counts)+geom_point(aes(distance,contact.close/exp(true_log_mean_cclose-true_log_decay)), alpha=0.1)+
    geom_line(aes(distance,exp(true_log_decay)))+scale_x_log10()+scale_y_log10()
  #statistics
  dset_statistics(biases,counts)
  return(list(biases=biases, counts=counts))
}

fill_zeros = function(biases,counts) {
  newcounts=CJ(biases[,paste(id,pos)],biases[,paste(id,pos)])
  newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
  newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
  newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
  newcounts=newcounts[pos1<pos2]
  setkey(newcounts, id1, id2, pos1, pos2)
  setkey(counts, id1, id2, pos1, pos2)
  newcounts=counts[newcounts]
  newcounts[is.na(contact.close),contact.close:=0]
  newcounts[is.na(contact.far),contact.far:=0]
  newcounts[is.na(contact.up),contact.up:=0]
  newcounts[is.na(contact.down),contact.down:=0]
  return(newcounts)
}




dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = .1, eRJ = 5, eDE = 7)
counts=dset$counts
biases=dset$biases

dset_statistics(biases,counts)
meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 10)

#### optimization wihout prior guesses
smfit = stan_model(file = "sparse_cs_norm_fit_meanfield.stan")
maxcount=3
system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, Krow=100, Kdiag=10,
                                      lambda_nu=1, lambda_delta=1, lambda_diag=.1, verbose = T))
#save(op, file = "data/caulo_all_op_lambda1.RData")
#compare offsets
c(op$par$eC,op$par$eRJ,op$par$eDE)
#compare decays
mcounts=melt(counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
             variable.name = "category", value.name = "count")[count>maxcount]
mcounts[,distance:=abs(pos2-pos1)]
mcounts=rbind(mcounts[,.(distance,count,category)],meanfield$Nkd[count<=maxcount,.(distance=mdist,count,category="meanfield")])
mcounts[category=="contact.close",fij:=exp(op$par$log_decay_close)]
mcounts[category=="contact.far",fij:=exp(op$par$log_decay_far)]
mcounts[category=="contact.up",fij:=exp(op$par$log_decay_up)]
mcounts[category=="contact.down",fij:=exp(op$par$log_decay_down)]
mcounts[category=="meanfield",fij:=exp(op$par$log_decay_mf)]
ggplot(mcounts[,.SD[sample(.N,min(.N,10000))],by=category])+scale_y_log10()+scale_x_log10()+
  geom_line(aes(distance, fij, colour=category))+geom_line(data=counts, aes(distance,exp(true_log_decay)))
#biases
#pbegin=35100000 
#pend=35200000 
#pbegin=3166716
#pend=3191637
#pbegin=73780165
#pend=74230165
pbegin=1
pend=10000
#nu
a=cbind(biases,data.table(opt=exp(op$par$log_nu)))
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op$par$eDE)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op$par$eDE)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op$par$eRJ)),colour="red")+
  geom_point(aes(pos, opt),colour="blue", shape=0)+
  geom_line(aes(pos, exp(true_log_nu)))
#delta relative to nu
a=cbind(biases,data.table(opt=exp(op$par$log_delta),
                          log_nu=op$par$log_nu))
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/(exp(log_nu+op$par$eDE))),colour="orange")+
  geom_point(aes(pos, dangling.R/(exp(log_nu+op$par$eDE))),colour="pink")+
  geom_point(aes(pos, rejoined/(exp(log_nu+op$par$eRJ))),colour="red")+
  geom_point(aes(pos, opt),colour="orange", shape=0)+
  geom_point(aes(pos, 1/opt),colour="pink", shape=0)+
  geom_line(aes(pos, exp(true_log_delta)))
#all 3 types
a=cbind(biases,data.table(RJ=exp(op$par$eRJ+op$par$log_nu),
                          DL=exp(op$par$eDE+op$par$log_nu+op$par$log_delta),
                          DR=exp(op$par$eDE+op$par$log_nu-op$par$log_delta)))
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, rejoined), colour="red")+ geom_point(aes(pos, RJ), shape=0, colour="red")+
  geom_point(aes(pos, dangling.L), colour="orange")+ geom_point(aes(pos, DL), shape=0, colour="orange")+
  geom_point(aes(pos, dangling.R), colour="pink")+ geom_point(aes(pos, DR), shape=0, colour="pink")+
  geom_line(aes(pos, exp(true_log_mean_DR)), colour="pink")+
  geom_line(aes(pos, exp(true_log_mean_DL)), colour="orange")+
  geom_line(aes(pos, exp(true_log_mean_RJ)), colour="red")

## predict on full dataset
smpred = stan_model(file = "sparse_cs_norm_predict.stan")
system.time(op.pred <- predict_full(smpred, biases, counts, op, Kdiag=10, verbose = T))
counts[,log_decay:=op.pred$par$log_decay]
counts[,log_mean_cup:=op.pred$par$log_mean_cup]
counts[,log_mean_cdown:=op.pred$par$log_mean_cdown]
counts[,log_mean_cfar:=op.pred$par$log_mean_cfar]
counts[,log_mean_cclose:=op.pred$par$log_mean_cclose]
#diagonal decay
ggplot(counts[sample(.N,min(.N,10000))])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(abs(pos2-pos1), contact.close/exp(log_mean_cclose-log_decay)), alpha=0.1) +
  geom_line(aes(abs(pos2-pos1), exp(log_decay)))
#residuals
ggplot(counts[abs(pos2-pos1)>1e4][sample(.N,min(.N,10000))])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(abs(pos2-pos1), contact.far/exp(log_mean_cfar)), alpha=0.1)


#correlations
biases[,log_nu:=op$par$log_nu]
biases[,log_delta:=op$par$log_delta]
biases[,log_mean_DL:=op$par$log_mean_DL]
biases[,log_mean_DR:=op$par$log_mean_DR]
#decay
ggplot(melt(counts[,.(distance,true_log_decay,log_decay)], id.vars = "distance"), aes(distance,value,colour=variable))+
  geom_line() + geom_point(size=0.5)
ggplot(counts)+geom_point(aes(true_log_decay,log_decay))+stat_function(fun=identity)
#cclose
ggplot(melt(counts[,.(.I,distance,true_log_mean_cclose,log_mean_cclose)], id.vars = c("I","distance")), aes(I,value,colour=variable))+
  geom_line() + geom_point(size=0.5) + xlim(0,100)
ggplot(counts)+geom_point(aes(true_log_mean_cclose,log_mean_cclose))+stat_function(fun=identity)
#nu
ggplot(melt(biases[,.(pos,true_log_nu,log_nu)], id.vars = "pos"), aes(pos,value,colour=variable))+
  geom_line() + geom_point(size=0.5) #+ xlim(1000,50000)
ggplot(biases)+geom_point(aes(true_log_nu,log_nu))+stat_function(fun=identity)
#delta
ggplot(melt(biases[,.(pos,true_log_delta,log_delta)], id.vars = "pos"), aes(pos,value,colour=variable))+
  geom_line() + geom_point(size=0.5) #+ xlim(1000,50000)
ggplot(biases)+geom_point(aes(true_log_delta,log_delta))+stat_function(fun=identity)
#dangling L
ggplot(melt(biases[,.(pos,true_log_mean_DL,log_mean_DL)], id.vars = "pos"), aes(pos,value,colour=variable))+
  geom_line() + geom_point(size=0.5) + xlim(1000,50000)
ggplot(biases)+geom_point(aes(true_log_mean_DL,log_mean_DL))+stat_function(fun=identity)
#dangling R
ggplot(melt(biases[,.(pos,true_log_mean_DR,log_mean_DR)], id.vars = "pos"), aes(pos,value,colour=variable))+
  geom_line() + geom_point(size=0.5) + xlim(1000,50000)
ggplot(biases)+geom_point(aes(true_log_mean_DR,log_mean_DR))+stat_function(fun=identity)



### effect of mean field approximation
biases=fread("data/rao_HICall_chr20_35000000-36000000_biases.dat")[id<1000] #1000
setkey(biases,id)
counts=fread("data/rao_HICall_chr20_35000000-36000000_counts.dat")[id1<1000&id2<1000]
counts=fill_zeros(biases,counts)

biases=fread("data/rao_HIC035_chr20_35000000-36000000_biases.dat") #212
setkey(biases,id)
counts=fread("data/rao_HIC035_chr20_35000000-36000000_counts.dat")
counts=fill_zeros(biases,counts)

biases=fread("data/caulo_3000000-4000000_biases.dat") #525
setkey(biases,id)
counts=fread("data/caulo_3000000-4000000_counts.dat")
counts=fill_zeros(biases,counts)

dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = .1, eRJ = 5, eDE = 7) #high coverage
dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = -2, eRJ = 1, eDE = 3) #low coverage

counts=dset$counts
biases=dset$biases
dset_statistics(biases, counts)
meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)
smfit = stan_model(file = "sparse_cs_norm_fit_meanfield.stan")
ops = foreach (maxcount=-1:10) %:% foreach (repetition=1:10) %dopar% {
  message("***** ",maxcount)
  #fit
  a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, Krow=100, Kdiag=10,
                                       lambda_nu=.3, lambda_delta=1.2, lambda_diag=1.6, verbose = T, iter=100000))
  op$par$time=a[1]+a[4]
  op$par$maxcount=maxcount
  op$par$repetition=repetition
  #predict
  op.pred <- predict_full(smpred, biases, counts, op, Kdiag=10, verbose = T)
  op$pred=op.pred$par
  op#ops[[paste0("op_",maxcount,"_",repetition)]]=op
}
ops = unlist(ops, recursive=F)
names(ops) <- seq(length(ops))
#save(ops, file = "data/rao_HICall_chr20_300k_mf.RData")
#load("data/caulo_mf_3-4M.RData")

#precision and accuracy on single-valued params
single_params = data.table(
  sapply(c("eC","eRJ","eDE","lambda_nu","lambda_delta","lambda_diag","alpha","deviance_proportion_explained",
           "repetition","maxcount","time"),
       function(y){sapply(names(ops), function(x){ops[[x]]$par[[y]]})}))
single_params[,maxcount:=as.character(maxcount)]
single_params[maxcount=="-1",maxcount:="exact"]
single_params[,maxcount:=ordered(maxcount,levels=unique(maxcount))]
ggplot(melt(single_params,id.vars=c("maxcount","repetition")))+geom_boxplot(aes(maxcount,value,colour=(maxcount=="exact")))+
  facet_wrap(~variable, scales="free_y")+guides(colour=F)+
  labs(title="Rao HICall chr20 300k stretch, 90% zeros", x="mean field threshold", y=NULL)
ggsave(filename = "images/rao_HICall_chr20_300k_mf_singleparams.png", width=10, height=7.5)

#precision on vector-valued params
all_params = data.table(
  sapply(names(ops[["1"]]$par), function(y){lapply(names(ops), function(x){ops[[x]]$par[[y]]})}))
all_params[,maxcount:=as.character(maxcount)]
all_params[maxcount=="-1",maxcount:="exact"]
all_params[,maxcount:=ordered(maxcount,levels=unique(maxcount))]
all_params[,repetition:=as.integer(repetition)]
mult_params=melt(all_params[,.(maxcount,repetition,log_nu,log_delta,beta_diag)], id.vars=c("maxcount","repetition"))
mult_params=mult_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("maxcount","repetition","variable")][
                         ,.(std=sd(value)),by=c("maxcount","idx","variable")]
ggplot(mult_params)+
  geom_boxplot(aes(maxcount,std,colour=(maxcount=="exact")))+facet_wrap(~variable)+
  guides(colour=F)+scale_y_log10()+labs(title="Rao HICall chr20 300k stretch, 90% zeros: precision", x="mean field threshold",
                                        y="standard deviation across repeats, for each coefficient")
ggsave(filename = "images/rao_HICall_chr20_300k_mf_multiparams_precision.png", width=10, height=7.5)

#accuracy on vector-valued params
ref_params=melt(all_params[maxcount=="exact",.(repetition,log_nu,log_delta,beta_diag)], id.vars="repetition")
ref_params=ref_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","repetition")][
  ,.(med=median(value)),keyby=c("idx","variable")]
mult_params=melt(all_params[,.(maxcount,repetition,log_nu,log_delta,beta_diag)], id.vars=c("maxcount","repetition"))
mult_params=mult_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("maxcount","repetition","variable")]
setkey(mult_params, idx, variable)
mult_params=ref_params[mult_params][,.(prec=dist(rbind(value,med))),by=c("maxcount","idx","variable")]
ggplot(mult_params)+
  geom_boxplot(aes(maxcount,prec,colour=(maxcount=="exact")))+facet_wrap(~variable)+
  guides(colour=F)+labs(title="Rao HICall chr20 300k stretch, 90% zeros: accuracy", x="mean field threshold",
                                        y="distance from median of exact calculation")+ylim(0,5)
ggsave(filename = "images/rao_HICall_chr20_300k_mf_multiparams_accuracy.png", width=10, height=7.5)

#plot nu
#pbegin=35100000 
#pend=35200000 
#pbegin=3166716
#pend=3266716
pbegin=0
pend=20000
mean_params=melt(all_params[,.(maxcount,repetition,log_nu,log_delta,eRJ,eDE)], id.vars=c("maxcount","repetition"))
mean_params=mean_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","maxcount","repetition")][
  ,.(med=median(value)),keyby=c("maxcount","variable","idx")]
setkey(mean_params, idx)
ggplot(biases[pos>=pbegin&pos<=pend])+scale_y_log10()+coord_cartesian(ylim=c(0.1,10))+
  geom_point(aes(pos, dangling.L/exp(mean_params[variable=='eDE'&maxcount=="exact",med])),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(mean_params[variable=='eDE'&maxcount=="exact",med])),colour="pink")+
  geom_point(aes(pos, rejoined/exp(mean_params[variable=='eRJ'&maxcount=="exact",med])),colour="red")+
  geom_line(data=biases[mean_params[variable=='log_nu']][pos>=pbegin&pos<=pend],aes(pos, exp(med), colour=maxcount))+
  geom_line(aes(pos, exp(true_log_nu)), colour="black")+
  #labs(title="Rao HICall chr20 300k stretch: 100k example for nu", x="genomic position", y="biases")
  labs(title="Caulobacter: 100k example for nu", x="genomic position", y="biases")
#ggsave(filename = "images/rao_HICall_chr20_300k_mf_nu_example.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-4M_mf_nu_example.png", width=10, height=7.5)

#plot delta
ggplot(biases[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(mean_params[variable=='eDE'&maxcount=="exact",med])),colour="orange")+
  geom_line(data=dcast(biases[mean_params], ...~variable, value.var="med")[pos>=pbegin&pos<=pend],
            aes(pos, exp(log_nu+log_delta), colour=maxcount))+
  #labs(title="Rao HICall chr20 300k stretch: 100k example for delta", x="genomic position", y="biases")
  labs(title="Caulobacter: 100k example for delta", x="genomic position", y="biases")
#ggsave(filename = "images/rao_HICall_chr20_300k_mf_delta_example.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-4M_mf_delta_example.png", width=10, height=7.5)

#plot decay
decay_params = cbind(
  data.table(sapply(c("maxcount","repetition"), function(y){sapply(names(ops), function(x){ops[[x]]$par[[y]]})})),
  data.table(sapply(names(ops[["1"]]$pred), function(y){lapply(names(ops), function(x){ops[[x]]$pred[[y]][1:counts[id1==1,.N]]})})))
decay_params[,maxcount:=as.character(maxcount)]
decay_params[maxcount=="-1",maxcount:="exact"]
decay_params[,maxcount:=ordered(maxcount,levels=unique(maxcount))]
decay_params[,repetition:=as.integer(repetition)]
decay_params=melt(decay_params, id.vars=c("maxcount","repetition"))
decay_params=decay_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","repetition","maxcount")]
decay_params=decay_params[,.(med=median(value)),keyby=c("idx","variable","maxcount")]
setkey(decay_params, idx)
decay_params[,dist:=rep(counts[id1==1,pos2-pos1],each=decay_params[idx==1,.N])]
decay_params=dcast(decay_params, ...~variable, value.var="med")
ggplot(decay_params)+scale_y_log10()+scale_x_log10()+
  geom_line(aes(abs(dist), exp(log_decay), colour=maxcount))+
  geom_line(data=counts[,.(decay=exp(true_log_decay)),by=distance],aes(distance,decay),colour="black")+
  #labs(title="Rao HICall chr20 300k stretch: diagonal decay", x="distance", y="decay")
  #labs(title="Caulobacter: 100k example for decay", x="genomic position", y="biases")
  labs(title="Toy dataset, high coverage, 40% zeros: decay", x="genomic position", y="biases")
#ggsave(filename = "images/rao_HICall_chr20_300k_mf_decay.png", width=10, height=7.5)
ggsave(filename = "images/toy_high_mf_decay_twoalphas_100.png", width=10, height=7.5)
#





### effect of removing zeros and ones
dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = .1, eRJ = 5, eDE = 7) #high coverage
dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = -2, eRJ = 1, eDE = 3) #low coverage
dset_statistics(dset$biases, dset$counts)
counts=dset$counts
biases=dset$biases
meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 50)
smfit = stan_model(file = "sparse_cs_norm_fit_meanfield.stan")
ops=list()
for (mincount in -1:2) {
  for (repetition in 1:10) {
    message("***** ",mincount)
    #fit
    a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=-1, Krow=100, Kdiag=10,
                                               lambda_nu=.3, lambda_delta=1.2, lambda_diag=1.6, verbose = T,
                                               mincount=mincount))
    op$par$time=a[1]+a[4]
    op$par$mincount=mincount
    op$par$repetition=repetition
    #predict
    op.pred <- predict_full(smpred, biases, counts, op, Kdiag=10, verbose = T)
    op$pred=op.pred$par
    ops[[paste0("op_",mincount,"_",repetition)]]=op
  }
}
save(ops, file = "data/toy_zero_high.RData")
#load("data/toy_mf_high.RData")
#precision and accuracy on single-valued params
single_params = data.table(
  sapply(c("eC","eRJ","eDE","lambda_nu","lambda_delta","lambda_diag","alpha","deviance_proportion_explained",
           "repetition","mincount","time"),
         function(y){sapply(names(ops), function(x){ops[[x]]$par[[y]]})}))
single_params[,mincount:=as.character(mincount)]
single_params[mincount=="-1",mincount:="exact"]
single_params[,mincount:=ordered(mincount,levels=unique(mincount))]
ggplot(melt(single_params[alpha<100000&lambda_diag<1000],id.vars=c("mincount","repetition")))+geom_boxplot(aes(mincount,value,colour=(mincount=="exact")))+
  facet_wrap(~variable, scales="free_y")+guides(colour=F)+
  labs(title="Toy dataset, high coverage, 75% zeros", x="mean field threshold", y=NULL)
ggsave(filename = "toy_zero_singleparams_low.png", width=10, height=7.5)
#precision on vector-valued params
all_params = data.table(
  sapply(names(ops$`op_-1_1`$par), function(y){lapply(names(ops), function(x){ops[[x]]$par[[y]]})}))
all_params[,mincount:=as.character(mincount)]
all_params[mincount=="-1",mincount:="exact"]
all_params[,mincount:=ordered(mincount,levels=unique(mincount))]
all_params[,repetition:=as.integer(repetition)]
mult_params=melt(all_params[,.(mincount,repetition,log_nu,log_delta,beta_diag)], id.vars=c("mincount","repetition"))
mult_params=mult_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("mincount","repetition","variable")][
  ,.(std=sd(value)),by=c("mincount","idx","variable")]
ggplot(mult_params)+
  geom_boxplot(aes(mincount,std,colour=(mincount=="exact")))+facet_wrap(~variable)+
  guides(colour=F)+scale_y_log10()+labs(title="Toy dataset, high coverage, 40% zeros: precision", x="mean field threshold",
                                        y="standard deviation across repeats, for each coefficient")
ggsave(filename = "toy_zero_multiparams_precision_high.png", width=10, height=7.5)
#accuracy on vector-valued params
ref_params=melt(all_params[mincount=="exact",.(repetition,log_nu,log_delta,beta_diag)], id.vars="repetition")
ref_params=ref_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","repetition")][
  ,.(med=median(value)),keyby=c("idx","variable")]
mult_params=melt(all_params[,.(mincount,repetition,log_nu,log_delta,beta_diag)], id.vars=c("mincount","repetition"))
mult_params=mult_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("mincount","repetition","variable")]
setkey(mult_params, idx, variable)
mult_params=ref_params[mult_params][,.(prec=dist(rbind(value,med))),by=c("mincount","idx","variable")]
ggplot(mult_params)+
  geom_boxplot(aes(mincount,prec,colour=(mincount=="exact")))+facet_wrap(~variable)+
  guides(colour=F)+labs(title="Toy dataset, high coverage, 40% zeros: accuracy", x="mean field threshold",
                        y="distance from median of exact calculation")
ggsave(filename = "toy_zero_multiparams_accuracy_high.png", width=10, height=7.5)

