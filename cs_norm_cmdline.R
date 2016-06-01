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
registerDoParallel(cores=30)

setwd("/home/yannick/simulations/cs_norm")

optimize_all_meanfield = function(model, biases, counts, meanfield, maxcount, bf_per_kb=1, bf_per_decade=5,
                                  iter=10000, verbose=T, init=0, ...) {
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  cclose=counts[contact.close>maxcount,.(id1,id2,distance,count=contact.close)]
  cfar=counts[contact.far>maxcount,.(id1,id2,distance,count=contact.far)]
  cup=counts[contact.up>maxcount,.(id1,id2,distance,count=contact.up)]
  cdown=counts[contact.down>maxcount,.(id1,id2,distance,count=contact.down)]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               dmin=dmin, dmax=dmax, Kdiag=Kdiag,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               Nl=mf$Nkl[,.N], Nkl_count=mf$Nkl[,count], Nkl_cidx=mf$Nkl[,id], Nkl_N=mf$Nkl[,N], Nkl_levels=mf$Nkl[,sum(diff(N)!=0)+1],
               Nr=mf$Nkr[,.N], Nkr_count=mf$Nkr[,count], Nkr_cidx=mf$Nkr[,id], Nkr_N=mf$Nkr[,N], Nkr_levels=mf$Nkr[,sum(diff(N)!=0)+1],
               Nd=mf$Nkd[,.N], Nkd_count=mf$Nkd[,count], Nkd_d=mf$Nkd[,mdist], Nkd_N=mf$Nkd[,N], Nkd_levels=mf$Nkd[,sum(diff(N)!=0)+1])
  if (data$Nl==0) data$Nkl_levels=0
  if (data$Nr==0) data$Nkr_levels=0
  if (data$Nd==0) data$Nkd_levels=0
  message("Mean field optimization")
  message("Krow        : ", Krow)
  message("Kdiag       : ", Kdiag)
  message("Biases      : ", biases[,.N])
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  message("Left counts : ", mf$Nkl[,.N], " (", data$Nkl_levels, " levels)")
  message("Right counts: ", mf$Nkr[,.N], " (", data$Nkr_levels, " levels)")
  message("Decay counts: ", mf$Nkd[,.N], " (", data$Nkd_levels, " levels)")
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

predict_all_meanfield = function(model, biases, counts, meanfield, opt, bf_per_decade=5, verbose=T) {
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag_centered=opt$par$beta_diag_centered, alpha=opt$par$alpha)
  message("Mean field prediction")
  message("Kdiag       : ", Kdiag)
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
}







args=commandArgs(trailingOnly = T)
prefix=args[1]
if (length(args)>=2) suffix=args[2] else suffix=NULL

message("normalization on ",prefix)

biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))
#meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)
#save(meanfield, file = paste0("data/",prefix,"_meanfield_100.RData"))
load(paste0("data/",prefix,"_meanfield_100.RData"))

#### optimization wihout prior guesses
smfit = stan_model(file = "cs_norm_fit.stan")
smpred = stan_model(file = "cs_norm_predict.stan")
maxcount=4
a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                         bf_per_decade=5, verbose = T, iter=100000)) #, tol_rel_grad=1e3, tol_rel_obj=1e3
op$pred=predict_all_meanfield(smpred, biases, counts, meanfield, op, verbose=T)$par
save(op, file=paste0("data/",prefix,"_op_maxcount_",maxcount,suffix,".RData"))


