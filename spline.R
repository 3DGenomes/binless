library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)

setwd("/home/yannick/simulations/spline_stan")

#generate data according to a function
generate_data = function(fun1, fun2, x1min=0, x1max=1, x2min=0, x2max=1, npoints=100, sd=1) {
  x1=seq(from=x1min, to=x1max, length.out=npoints)
  x2=seq(from=x2min, to=x2max, length.out=npoints)
  f1=sapply(x1,fun1)
  f2=sapply(x2,fun2)
  #y=rnorm(npoints, mean=f1+f2, sd=sd)
  y=rnbinom(npoints, mu=f1+f2, size=sd)
  return(data.table(x1=x1,x2=x2,f1=f1,f2=f2,y=y))
}

stan_matrix_to_datatable = function(opt, x1, x2) {
  vals=data.table(opt)
  vals[,x1:=x1]
  vals[,x2:=x2]
  melt(data.table(vals), id.vars=c("x1","x2"))
}


#plot initial data
data=generate_data(fun1=function(x){1+10*exp(-(x-1.5)**2*2)}, fun2=function(x){5+5*sin(5*x)},
                   sd=5, x1min=0,x1max=3, x2min=-1, x2max=1, npoints=1000)
#data=data[y>0]
#data=generate_data(fun=function(x){5*sin(3*x)+5.1}, sd=1, xmin=0,xmax=3, npoints=500)
ggplot(data)+geom_point(aes(x1,y))+geom_line(aes(x1,f1))#+scale_y_log10()
ggplot(data)+geom_point(aes(x2,y))+geom_line(aes(x2,f2))#+scale_y_log10()
ggplot(data)+geom_point(aes(f1+f2,y))

#fit it with stan
sm = stan_model(file = "spline.stan")
op = optimizing(sm, data = list(N=data[,.N], K1=10, K2=10, y=data[,y], x1=data[,x1], x2=data[,x2]),
                as_vector=F, hessian=F, iter=10000)
#sf = stan(file="spline.stan", data = list(N=data[,.N], K=20, y=data[,y], x=data[,x]))
#launch_shinystan(sf)

c(op$par$lambda,op$par$edf)
#mat=stan_matrix_to_datatable(op$par$designmat, data[,x])
mat=stan_matrix_to_datatable(op$par$weighted, data[,x1], data[,x2])
#mat=stan_matrix_to_datatable(op$par$dofmat, 1:10)
data[,pred:=op$par$pred]
data[,s1:=op$par$splines[1,]]
data[,s2:=op$par$splines[2,]]
#compare to gam
fit=gam(data=data, formula = y~s(x1)+s(x2), family=gaussian())
#fit=gam(data=data, formula = y~s(x,bs="ps", m=c(2,2), k=10), family=nb())
data[,gam:=fit$fitted.values]
data[,gam1:=predict(fit, data, type="terms", terms="s(x1)")]
data[,gam2:=predict(fit, data, type="terms", terms="s(x2)")]
#plot result
ggplot(data)+geom_point(aes(x1,y),alpha=0.2)+geom_line(aes(x1,f1+f2),colour="black")+#ylim(0,40)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(x1,gam),colour="blue")+
  geom_line(aes(x1,pred),colour="red")#+scale_y_log10()
ggplot(data)+geom_point(aes(x2,y),alpha=0.2)+geom_line(aes(x2,f1+f2),colour="black")+#ylim(0,40)+
  #geom_line(data=mat, aes(x2,value,colour=variable))+
  geom_line(aes(x2,gam),colour="blue")+
  geom_line(aes(x2,pred),colour="red")#+scale_y_log10()
ggplot(data)+geom_point(aes(x1,y-f2),alpha=0.2)+geom_line(aes(x1,f1),colour="black")+#ylim(0,40)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(x1,gam1),colour="blue")+
  geom_line(aes(x1,s1),colour="red")#+scale_y_log10()
ggplot(data)+geom_point(aes(x2,y-f1),alpha=0.2)+geom_line(aes(x2,f2),colour="black")+#ylim(0,40)+
  #geom_line(data=mat, aes(x1,value,colour=variable))+
  geom_line(aes(x2,gam2),colour="blue")+
  geom_line(aes(x2,s2),colour="red")#+scale_y_log10()


