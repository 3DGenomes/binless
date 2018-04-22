compute_means_ori = function(cs, counts) {
  #compute background
  init=cs@par
  cpos=copy(counts)
  bsub=cs@biases[,.(id)] 
  bsub[,c("log_iota","log_rho"):=list(init$log_iota,init$log_rho)]
  cpos=merge(bsub[,.(id1=id,log_iota,log_rho)],cpos,by="id1",all.x=F,all.y=T)
  cpos=merge(bsub[,.(id2=id,log_iota,log_rho)],cpos,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
  cpos=merge(cbind(cs@design[,.(name)],eC=init$eC), cpos, by="name",all.x=F,all.y=T)
  cpos[,c("bin1","bin2","dbin"):=
         list(cut(pos1, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(pos2, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cpos=merge(cpos,init$decay[,.(name,dbin,log_decay)],by=c("name","dbin"))
  #compute signal
  if (cs@par$signal[,.N]>0 && length(cs@settings$sbins)>2) {
    signal = binless:::get_signal_matrix(cs, resolution = cs@settings$base.res, groups=cs@experiments[,.(name,groupname=name)])
    signal = rbind(signal[,.(name,bin1,bin2,phi)],signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
    cpos = signal[cpos,,on=c("name","bin1","bin2")]
  } else {
    cpos[,phi:=0]
  } 
  #assemble
  cpos[,log_mu.base:=eC + log_decay + phi]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  
  cpos=cpos[,.(id1,id2,name,pos1,pos2,distance,contact.close,contact.far,contact.up,contact.down,
               log_decay,lmu.close,lmu.far,lmu.up,lmu.down)]
  setkeyv(cpos,key(cs@counts))
  cpos
}   

gauss_dispersion_ori = function(cs, counts, weight=cs@design[,.(name,wt=1)], verbose=T) {
  if (verbose==T) cat(" Dispersion\n")
  #predict all means and put into table
  counts = compute_means_ori(cs,counts)
  stopifnot(cs@biases[,.N]==length(cs@par$log_iota))
  #
  #fit dispersion and exposures
  if (verbose==T) cat("  predict\n")
  bbegin=c(1,cs@biases[,.(name,row=.I)][name!=shift(name),row],cs@biases[,.N+1])
  cbegin=c(1,counts[,.(name,row=.I)][name!=shift(name),row],counts[,.N+1])
  data = list( Dsets=cs@design[,.N], Biases=cs@design[,uniqueN(genomic)], Decays=cs@design[,uniqueN(decay)],
               XB=as.array(cs@design[,genomic]), XD=as.array(cs@design[,decay]),
               SD=cs@biases[,.N], bbegin=bbegin,
               cutsitesD=cs@biases[,pos], rejoined=cs@biases[,rejoined],
               danglingL=cs@biases[,dangling.L], danglingR=cs@biases[,dangling.R],
               N=counts[,.N], cbegin=cbegin,
               counts_close=counts[,contact.close], counts_far=counts[,contact.far],
               counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               weight=as.array(weight[,wt]),
               log_iota=cs@par$log_iota, log_rho=cs@par$log_rho,
               log_mean_cclose=counts[,lmu.close], log_mean_cfar=counts[,lmu.far],
               log_mean_cup=counts[,lmu.up], log_mean_cdown=counts[,lmu.down])
  init=list(eC_sup=as.array(counts[,log(mean(contact.close/exp(lmu.close))),by=name][,V1]),
            eRJ=as.array(cs@biases[,.(name,frac=rejoined/exp((cs@par$log_iota+cs@par$log_rho)/2))][,log(mean(frac)),by=name][,V1]),
            eDE=as.array(cs@biases[,.(name,frac=(dangling.L/exp(cs@par$log_iota)+dangling.R/exp(cs@par$log_rho))/2)][
              ,log(mean(frac)),by=name][,V1]))
  init$mu=mean(exp(init$eC_sup[1]+counts[name==name[1],lmu.close]))
  init$alpha=max(0.001,1/(var(counts[name==name[1],contact.close]/init$mu)-1/init$mu))
  init$mu=NULL
  out=capture.output(op<-optimize_stan_model(model=binless:::stanmodels$gauss_dispersion, tol_param=cs@par$tol_disp,
                                             data=data, iter=cs@settings$iter, verbose=verbose, init=init,
                                             init_alpha=1e-9))
  #restrict tolerance if needed
  precision = max(abs(c(op$par["alpha"],recursive=T) - c(cs@par["alpha"],recursive=T)))
  cs@par$tol_disp = min(cs@par$tol_disp, max(cs@settings$tol, precision/10))
  #update parameters
  cs@par=modifyList(cs@par, op$par["alpha"])
  #cs@par$eC=cs@par$eC+op$par$eC_sup
  #
  #compute log-posterior
  Krow=cs@settings$Krow
  Kdiag=cs@settings$Kdiag
  cs@par$value = op$value + (Krow-2)/2*sum(log(cs@par$lambda_iota/exp(1))+log(cs@par$lambda_rho/exp(1))) +
    (Kdiag-2)/2*sum(log(cs@par$lambda_diag/exp(1)))
  if (verbose==T) {
    cat("  fit: dispersion",cs@par$alpha,"\n")
    cat("  log-likelihood = ",cs@par$value,"\n")
  }
  return(cs)
}

cts.common = binless:::gauss_common_muhat_mean(cs.new.0, cs.new.0@zeros, cs.new.0@settings$sbins)
biasmat = binless:::gauss_genomic_muhat_mean(cs.new.0, cts.common)
a = gauss_genomic_muhat_mean_ori(cs.ori.0, cts.common)
biasmat.ori=rbindlist(a)
biasmats=rbindlist(list(ori=biasmat.ori[,.(pos,cat,etahat,std)],new=biasmat[,.(pos,cat,etahat,std)]),use=T,id="origin")
dcast(biasmats,pos+cat~origin,value.var = "etahat")[,summary(new-ori)]
dcast(biasmats,pos+cat~origin,value.var = "std")[,summary(new-ori)]

op.new = binless:::gauss_genomic_optimize(biasmat, cs.new.0@design, cs.new.0@settings$Krow, cs.new.0@settings$sbins,
                                      cs.new.0@par$lambda_iota, cs.new.0@par$lambda_rho, verbose=T,
                                      max_perf_iteration=cs.new.0@settings$iter,
                                      convergence_epsilon=cs.new.0@par$tol_genomic,
                                      constrain=T)

generate_spline_base = binless:::generate_spline_base
op.ori = gauss_genomic_optimize_ori(a$bts, a$cts, cs.ori.0@biases, cs.ori.0@design, cs.ori.0@settings$Krow, cs.ori.0@settings$sbins,
                                          cs.ori.0@par$lambda_iota, cs.ori.0@par$lambda_rho, verbose=T,
                                          max_perf_iteration=cs.ori.0@settings$iter,
                                          convergence_epsilon=cs.ori.0@par$tol_genomic,
                                          constrain=T)

biasmats=rbindlist(list(ori=op.ori$biases[,.(pos,cat,etahat,std,eta)],new=op.new$biases[,.(pos,cat,etahat,std,eta)]),use=T,id="origin")
dcast(biasmats,pos+cat~origin,value.var = "etahat")[,summary(new-ori)]
dcast(biasmats,pos+cat~origin,value.var = "std")[,summary(new-ori)]
dcast(biasmats,pos+cat~origin,value.var = "eta")[,summary(new-ori)]
ggplot(biasmats)+geom_line(aes(pos,eta,colour=origin))+facet_wrap(~cat)
