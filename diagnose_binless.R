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
binless:::has_converged(cs)

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
ggplot(signals[name==unique(name)[2]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=beta))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
ggplot(signals[step>=step[.N]-1])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_grid(step~ name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()
ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=beta))+geom_raster(aes(bin2,bin1,fill=pmin(phihat,max(beta))))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()

#signal differential
ggplot(signals[name==unique(name)[1],.(step,value=beta-shift(beta)),by=c("name","bin1","bin2")])+geom_raster(aes(bin1,bin2,fill=value))+facet_wrap(~step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")+coord_fixed()

#plot last signal row (virtual 4C)
ggplot(signals[bin1==min(bin1)])+geom_line(aes(bin2,phi,group=name))+geom_point(aes(bin2,phihat),alpha=0.1)+facet_grid(step~name)+ylim(-2.5,10)

#biases vs step
biases.ori=merge(rbindlist(lapply(cs.ori@diagnostics$params[leg=="bias",log_iota],function(x){cs@biases[,.(name,pos,log_iota=x)]}),use=T,id="step"),
             rbindlist(lapply(cs.ori@diagnostics$params[leg=="bias",log_rho],function(x){cs@biases[,.(name,pos,log_rho=x)]}),use=T,id="step"))[name==name[1]]
biases.new=rbindlist(cs.new@diagnostics$params[leg=="bias",biases],use=T,id="step")[group==group[1]&cat%in%c("contact L","contact R")]
biases.new=merge(biases.new[cat=="contact L",.(step,pos,log_iota=eta)],biases.new[cat=="contact R",.(step,group,pos,log_rho=eta)],by=c("step","pos"))
biases=melt(rbindlist(list(new=biases.new[,.(step,pos,log_iota,log_rho)],ori=biases.ori[,.(step,pos,log_iota,log_rho)]),use=T,idcol="origin"),id=c("origin","step","pos"))
ggplot(biases)+geom_line(aes(pos,value,colour=origin))+facet_grid(step~variable)

#biases differential
ggplot(biases[,.(step,value=value-shift(value)),by=c("name","pos","variable")])+geom_line(aes(pos,value,colour=variable))+facet_grid(step~variable)

#decay vs step (x log scale)
decays.ori=rbindlist(cs.ori@diagnostics$params[leg=="decay",decay],idcol = "step", use=T)
decays.ori[,group:=factor(unclass(name))]
decays.ori[,name:=NULL]
decays.new=rbindlist(cs.new@diagnostics$params[leg=="decay",decay],idcol = "step", use=T)
decays.new[,group:=factor(group)]
decays=rbindlist(list(ori=decays.ori,new=decays.new),use=T,id="origin")
ggplot(decays[,.(origin,distance,log_decay,std,kappahat),by=c("group","step")])+
  geom_line(aes(distance,log_decay,group=origin,colour=origin))+geom_point(aes(distance,kappahat,colour=origin),alpha=0.1)+
  geom_errorbar(aes(distance,ymin=kappahat-std,ymax=kappahat+std,colour=origin),alpha=0.1)+facet_grid(group~step)+scale_x_log10()

#decay vs step (x linear scale)
ggplot(decays[,.(distance,log_decay,std,kappahat),by=c("group","step")])+
  geom_line(aes(distance,log_decay,group=group))+geom_point(aes(distance,kappahat,colour=group),alpha=0.1)+
  geom_errorbar(aes(distance,ymin=kappahat-std,ymax=kappahat+std,colour=group),alpha=0.1)+facet_wrap(~step)

#fit last decay to polymer model
decay=cs@par$decay[name==name[1],.(ldist=log(distance),ldec=log_decay,distance,decay=exp(log_decay))]
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



