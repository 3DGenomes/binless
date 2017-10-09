library(csnorm)
library(ggplot2)
library(data.table)
library(scales)

setwd("/Users/yannick/Documents/simulations/cs_norm")

load("foxp1ext_observed.RData")

nouter=20
lam2=5
tol_val=1e-2
#out=csnorm:::fast_binless(mat, mat[,nlevels(bin1)], nouter, lam2, tol_val)
#save(out,file="out.dat")

load("out.dat")
diff=as.data.table(csnorm:::fast_binless_difference(out, lam2, tol_val, 1))
save(out,diff,file="out.dat")

#to debug:
#dsymutil /Library/Frameworks/R.framework/Versions/3.3/Resources/library/csnorm/libs/csnorm.so
#R -d "valgrind --dsymutil=yes" -f fast_binless.R 2>&1 | tee valgrind.out


if (F) {
  #TODO: smooth decay
  #TODO: difference calculation
  #TODO: negative binomial
  a=as.data.table(out$mat)
  #observed matrix
  ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=log(observed)))+geom_raster(aes(bin2,bin1,fill=log(observed)))+
    facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"), na.value = "white")
  #expected matrix
  ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=log(expected)))+geom_raster(aes(bin2,bin1,fill=log(expected)))+
    facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
  #fitted biases
  ggplot(data.table(bin=1:nlevels(mat[,bin1]),log_biases=out$log_biases))+geom_point(aes(bin,log_biases,colour="cpp"))#+
    geom_point(data=rbind(mat[,.(bin1,bin2,observed)],mat[bin2!=bin1,.(bin2=bin1,bin1=bin2,observed)]
                          )[observed>0,mean(log(observed)),by=bin1][,.(bin1,V1-mean(V1))],aes(bin1,V2,colour="true"))
  #biases matrix
  ggplot(out$mat)+geom_raster(aes(bin1,bin2,fill=biases))+geom_raster(aes(bin2,bin1,fill=biases))+
    facet_wrap(~name)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
  #fitted decay
  ggplot(data.table(distance=1:nlevels(mat[,bin1]),log_decay=out$log_decay))+geom_point(aes(distance,log_decay,colour="cpp"))#+
    geom_point(data=mat[observed>0,mean(log(observed)),by=bin2-bin1][,.(bin=bin2,log_decay=V1-mean(V1))],aes(bin,log_decay,colour="true"))
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
  
  #diagnostics
  b2=rbindlist(out$diagnostics,use=T,idcol = "step")
  ggplot(b[,.SD[1],by=c("step","bin1"),.SDcols="biases"][bin1>10&bin1<40])+geom_line(aes(step,biases,colour=bin1,group=bin1))
  ggplot(b[,.(step,d=bin2-bin1,decay)][,.SD[1],by=c("step","d")][d<100])+geom_line(aes(d,decay,colour=step,group=step))
  ggplot(b[name==name[1]])+geom_raster(aes(bin1,bin2,fill=signal))+geom_raster(aes(bin2,bin1,fill=signal))+
    facet_wrap(~step)+coord_fixed()+scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))
  
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
  
  
  iter=diff[,.(name,bin1,bin2,observed,I=.I-1)]
  a=fread("test.dat")
  a=dcast(a, V2~V1)
  setnames(a,"V2","I")
  iter=merge(a,iter,by="I")

  #phihat diff vs signal
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(phihat,5)),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(phihat,5)),data=out$mat)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  #phihat_ref diff vs phihat signal
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(phihat_ref,5)),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(phihat,5)),data=out$mat)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  #weights diff vs signal
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(weights,5)),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(weights,5)),data=out$mat)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  #weights_ref vs weights
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(weights,5)),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(weights_ref,5)),data=iter)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  #phihat vs deltahat
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(phihat,5)),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(deltahat,5)),data=iter)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  #delta vs deltahat
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=delta),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(deltahat,5)),data=iter)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  #phi_ref vs phihat_ref
  ggplot(diff)+geom_raster(aes(bin1,bin2,fill=pmin(phi_ref,5)),data=iter)+
    geom_raster(aes(bin2,bin1,fill=pmin(phihat_ref,5)),data=iter)+facet_wrap(~name)+
    scale_fill_gradient2(low=muted("blue"),mid="white",high=muted("red"))+coord_fixed()
  
}
