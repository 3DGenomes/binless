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

convert_to_simple = function(data) {
  #convert full cs norm to simple cs norm input
  newdata=data[,.(begin1,begin2,N,rejoined.1,rejoined.2,dangling.L.1,dangling.L.2,dangling.R.1,dangling.R.2,distance)]
  return(newdata[,.(N=sum(N)),by=c("begin1","begin2","rejoined.1","rejoined.2","dangling.L.1","dangling.L.2",
                                   "dangling.R.1","dangling.R.2","distance")])
}

#load data
data=fread("cs_norm_binned_rao_HiC035_chr22.dat", stringsAsFactors = T)
data=fread("cs_norm_binned_rao_HiC036_chr22.dat", stringsAsFactors = T)
data=fread("cs_norm_binned_rao_HiC036_chr22_10k.dat", stringsAsFactors = T)
data=fread("cs_norm_binned_rao_HiCall_chr22.dat", stringsAsFactors = T)
data=fread("cs_norm_binned_caulo.dat", stringsAsFactors = T)

data=fread("cs_norm_caulo.dat", stringsAsFactors = T)

data[,summary(N)]
ggplot(data)+geom_point(aes(distance,N),alpha=0.01)+
  scale_x_log10()+scale_y_log10()
data[,random.1:=random.1+1]
data[,random.2:=random.2+1]


#fit it with stan and gam
sm = stan_model(file = "binned_cs_norm.stan")
subdata=data[,.SD[sample(.N,10000)]]
system.time(subop <- optimizing(sm, data = list(N=subdata[,.N], K=10, B=4, #dangling L/R and rejoined
                                                count=subdata[,N], type=subdata[,as.integer(category)],
                                                dist=subdata[,distance],
                                                biases=as.matrix(subdata[,.(rejoined.1,rejoined.2,
                                                                      dangling.L.1,dangling.L.2,
                                                                      dangling.R.1,dangling.R.2,
                                                                      random.1, random.2)])),
                             as_vector=F, hessian=F, iter=10000))
system.time(op <- optimizing(sm, data = list(N=data[,.N], K=10, B=4, #dangling L/R and rejoined
                                             count=data[,N], type=data[,as.integer(category)],
                                             dist=data[,distance],
                                             biases=as.matrix(data[,.(rejoined.1,rejoined.2,
                                                                      dangling.L.1,dangling.L.2,
                                                                      dangling.R.1,dangling.R.2,
                                                                      random.1, random.2)])),
                             init = subop$par, 
                             as_vector=F, hessian=F, iter=10000))


#sf = stan(file="spline.stan", data = list(N=data[,.N], K=20, y=data[,y], x=data[,x]))
#launch_shinystan(sf)
#system.time(fit <- gam(N ~ s(log(distance)) + log(rejoined.1) + log(rejoined.2), data=data, family=nb()))
system.time(fit4 <- gam(N ~ s(log(distance),bs="ps", m=c(2,2), k=10) +
                         category:(log(dangling.L.1) + log(dangling.R.1) + log(rejoined.1) +
                         log(dangling.L.2) + log(dangling.R.2) + log(rejoined.2) + log(random.1) + log(random.2)),
                       data=data, family=nb()))
#system.time(fit <- scam(N ~ s(log(distance), bs="mpd", m=2, k=10) + log(rejoined.1) + log(rejoined.2),
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
ggplot(data)+geom_point(aes(distance,N/gam*gamfij),alpha=0.01)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(distance,fij),colour="red")+scale_x_log10()+scale_y_log10()+
  geom_line(aes(distance,gamfij),colour="blue")

ggplot(data)+geom_point(aes(distance,N/gam*gamfij),alpha=0.01)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  scale_x_log10()+scale_y_log10()+
  geom_line(aes(distance,gamfij),colour="blue")


ggplot(data)+geom_point(aes(distance,N),alpha=0.01)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(distance,pred),colour="red")+scale_x_log10()+
  geom_line(aes(distance,gam),colour="blue")

ggplot(data[sample(.N,10000)])+geom_point(aes(distance,N*fij/pred),alpha=0.2, colour="red")+scale_x_log10()+scale_y_log10()+
  geom_point(aes(distance,N*gamfij/gam),alpha=0.2, colour="blue")
  
ggplot(data[sample(.N,10000)])+geom_point(aes(fij,gamfij),alpha=0.2, colour="red")+scale_x_log10()+scale_y_log10()+
  stat_function(fun=function(x){x})




  