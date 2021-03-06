library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)

#This script is to group all commands related to diagnosing a run that does not converge

#check convergence and see precision
binless:::has_converged(cs, "signal")

#runtime: total and vs step
cs@diagnostics$params[,sum(runtime)/3600]
ggplot(cs@diagnostics$params[,.(step,leg,runtime)])+geom_line(aes(step,runtime,colour=leg))+scale_y_log10()

#log-likelihoods and parameters vs step
plot_diagnostics(cs)$plot
plot_diagnostics(cs)$plot2

#signal
plot_binless_matrix(cs@par$signal,upper="phi",lower="beta",trans="identity")

#signal vs step
signals=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  if ("signal" %in% cs@diagnostics$params[step==i,leg]) {
    sig=copy(cs@diagnostics$params[step==i&leg=="signal",signal][[1]])
    sig[,step:=i]
    sig
  }
}
ggplot(signals[name==unique(name)[1]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
ggplot(signals[name==unique(name)[1]])+geom_raster(aes(bin1,bin2,fill=log(weight)))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
ggplot(signals[name==unique(name)[2]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
ggplot(signals[step>=step[.N]-1])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_grid(step~ name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=beta))+geom_raster(aes(bin2,bin1,fill=pmin(phihat,max(beta))))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()

#signal differential
ggplot(signals[name==unique(name)[1],.(step,value=beta-shift(beta)),by=c("name","bin1","bin2")])+geom_raster(aes(bin1,bin2,fill=value))+facet_wrap(~step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()

#plot last signal row (virtual 4C)
ggplot(signals[bin1==min(bin1)])+geom_line(aes(bin2,phi,group=name))+geom_point(aes(bin2,phihat),alpha=0.1)+facet_grid(step~name)+ylim(-2.5,10)

#biases vs step
biases=rbindlist(cs@diagnostics$params[leg=="bias",biases],use=T,id="step")[cat%in%c("contact L","contact R")]
biases=rbind(biases[cat=="contact L",.(step,bias="log_iota",group,pos,eta,etahat,std,nobs)],
             biases[cat=="contact R",.(step,bias="log_rho",group,pos,eta,etahat,std,nobs)])
ggplot(biases)+geom_line(aes(pos,eta))+geom_point(aes(pos,etahat),alpha=0.1)+
  facet_grid(step~group+bias)

#biases differential
ggplot(biases[,.(step,eta=eta-shift(eta)),by=c("group","pos","bias")])+geom_line(aes(pos,eta,colour=bias))+facet_grid(step~group+bias)

#decay vs step (x log scale)
decays=rbindlist(cs@diagnostics$params[leg=="decay",decay],idcol = "step", use=T)
decays[,group:=factor(group)]
ggplot(decays)+
  geom_line(aes(distance,log_decay,group=group))+geom_point(aes(distance,kappahat,colour=group),alpha=0.1)+
  geom_errorbar(aes(distance,ymin=kappahat-std,ymax=kappahat+std,colour=group),alpha=0.1)+facet_grid(group~step)+scale_x_log10()

#decay vs step (x linear scale)
ggplot(decays)+
  geom_line(aes(distance,log_decay,group=group))+geom_point(aes(distance,kappahat,colour=group),alpha=0.1)+
  geom_errorbar(aes(distance,ymin=kappahat-std,ymax=kappahat+std,colour=group),alpha=0.1)+facet_grid(group~step)

#fit last decay to polymer model
decay=cs@par$decay[group==group[1],.(ldist=log(distance),ldec=log_decay,distance,decay=exp(log_decay))]
fit=lm(ldec~ldist,data=decay[distance>1e4&distance<5e5])
summary(fit)
decay[,model:=exp(fit$coefficients[["(Intercept)"]]+fit$coefficients[["ldist"]]*ldist)]
ggplot(decay)+geom_line(aes(distance,decay,colour="decay"))+
  scale_x_log10()+scale_y_log10()+geom_line(aes(distance,model,colour="model"))

#residuals at last step (along first row of the matrix)
a=cs@diagnostics$residuals[step==max(step)]
ggplot()+geom_point(aes(unclass(bin),count/nobs),alpha=0.5,data=a)+facet_wrap(~name)+scale_y_log10()+
  geom_line(aes(unclass(bin),value/nobs,colour=variable),
            data=melt(a[,.(name,bin,signal,decay,bias,mean,nobs)],id=c("name","bin","nobs")))



