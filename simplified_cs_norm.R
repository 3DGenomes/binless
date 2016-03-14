library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)
library(scam)

setwd("/home/yannick/simulations/spline_stan")

stan_matrix_to_datatable = function(opt, x1, x2) {
  vals=data.table(opt)
  vals[,x1:=x1]
  vals[,x2:=x2]
  melt(data.table(vals), id.vars=c("x1","x2"))
}


#load data
#data=fread("cs_norm_caulo.dat")
#data=fread("cs_norm_rao_HiC035_chr22.dat")
data=fread("cs_norm_binned_rao_HiC035_chr22.dat")
#data[,is.close:=factor(is.close)]

#fit it with stan and gam
sm = stan_model(file = "simplified_cs_norm.stan")
subdata=data[,.SD[sample(.N,10000)],by=is.close]
system.time(subop <- optimizing(sm, data = list(N=subdata[,.N], K=10, count=subdata[,N], dist=subdata[,distance],
                                             isclose=subdata[,as.integer(is.close=="TRUE")],
                                             rejoined1=subdata[,rejoined.1], rejoined2=subdata[,rejoined.2]),
                             as_vector=F, hessian=F, iter=10000))
system.time(op <- optimizing(sm, data = list(N=data[,.N], K=10, count=data[,N], dist=data[,distance],
                                             isclose=data[,as.integer(is.close=="TRUE")],
                                             rejoined1=data[,rejoined.1], rejoined2=data[,rejoined.2]),
                             init = subop$par,
                             as_vector=F, hessian=F, iter=10000))


#sf = stan(file="spline.stan", data = list(N=data[,.N], K=20, y=data[,y], x=data[,x]))
#launch_shinystan(sf)
#system.time(fit <- gam(N ~ s(log(distance)) + is.close + log(rejoined.1) + log(rejoined.2), data=data, family=nb()))
system.time(fit <- gam(N ~ s(log(distance),bs="ps", m=c(2,2), k=10) + is.close + log(rejoined.1) + log(rejoined.2),
                       data=data, family=nb()))
#system.time(fit <- scam(N ~ s(log(distance), bs="mpd", m=2, k=10) + is.close + log(rejoined.1) + log(rejoined.2),
#                       data=data, family=negbin(theta=1.93)))

c(op$par$lambda,op$par$edf)
#mat=stan_matrix_to_datatable(op$par$designmat, data[,x])
#mat=stan_matrix_to_datatable(op$par$weighted, data[,x1], data[,x2])
#mat=stan_matrix_to_datatable(op$par$dofmat, 1:10)
data[,pred:=op$par$pred]
data[,fij:=op$par$decay]
#compare to gam
#fit=gam(data=data, formula = y~s(x,bs="ps", m=c(2,2), k=10), family=nb())
data[,gam:=fit$fitted.values]
data[,gamfij:=exp(predict(fit, data, type="terms", terms="s(log(distance))"))]
#data[,gamfij:=exp(predict(fit, data, type="terms")[,"s(log(distance))"])]
#plot result
ggplot(data)+geom_point(aes(distance,N/pred*fij),alpha=0.01)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(distance,fij),colour="red")+scale_x_log10()+scale_y_log10()+
  geom_line(aes(distance,gamfij),colour="blue")
ggplot(data)+geom_point(aes(distance,N/gam*gamfij),alpha=0.05)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(distance,fij),colour="red")+scale_x_log10()+scale_y_log10()+
  geom_line(aes(distance,gamfij),colour="blue")

ggplot(data)+geom_point(aes(distance,N),alpha=0.01)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(distance,pred),colour="red")+scale_x_log10()+
  geom_line(aes(distance,gam),colour="blue")

ggplot(data[sample(.N,10000)])+geom_point(aes(distance,N*fij/pred),alpha=0.2, colour="red")+scale_x_log10()+scale_y_log10()+
  geom_point(aes(distance,N*gamfij/gam),alpha=0.2, colour="blue")
  
ggplot(data[sample(.N,10000)])+geom_point(aes(fij,gamfij),alpha=0.2, colour="red")+scale_x_log10()+scale_y_log10()+
  stat_function(fun=function(x){x})




  