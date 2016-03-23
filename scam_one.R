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

#generate data according to a function
generate_data = function(fun=sin, xmin=0, xmax=1, npoints=100, sd=1) {
  x=seq(from=xmin, to=xmax, length.out=npoints)
  f=sapply(x,fun)
  y=rnorm(npoints, mean=f, sd=sd)
  #y=rnbinom(npoints, mu=f, size=sd)
  return(data.table(x=x,f=f,y=y))
}

stan_matrix_to_datatable = function(opt, x) {
  vals=data.table(opt)
  vals[,x:=x]
  melt(data.table(vals), id.vars="x")
}


#plot initial data
data=generate_data(fun=function(x){1-6*x+0.5*x^3}, sd=3, xmin=0,xmax=3, npoints=1000)
#data=generate_data(fun=function(x){1-6*x}, sd=3, xmin=0,xmax=3, npoints=1000)
#data[,c("y","f"):=list(y-mean(f),f-mean(f))]
#data=data[y>0]
#data=generate_data(fun=function(x){5*sin(3*x)+5.1}, sd=1, xmin=0,xmax=3, npoints=500)
ggplot(data)+geom_point(aes(x,y))+geom_line(aes(x,f))#+scale_y_log10()

#fit it with stan
scm = stan_model(file = "scam_one_centered_params.stan")
scmi = stan_model(file = "scam_one_onlyspline.stan")
op=NULL
opi=NULL
op = optimizing(scm, data = list(N=data[,.N], K=10, y=data[,y], x=data[,x]),
                as_vector=F, hessian=F, iter=1000)
opi = optimizing(scmi, data = list(N=data[,.N], K=10, y=data[,y], x=data[,x]),
                as_vector=F, hessian=F, iter=1000)
#sf = stan(file="spline.stan", data = list(N=data[,.N], K=20, y=data[,y], x=data[,x]))
#launch_shinystan(sf)

c(op$par$lambda,op$par$edf,op$value)
c(opi$par$lambda,opi$par$edf,opi$value)
#mat=stan_matrix_to_datatable(op$par$designmat, data[,x])
mat=stan_matrix_to_datatable(opi$par$weighted, data[,x])
#mat=stan_matrix_to_datatable(op$par$dofmat, 1:10)
data[,pred:=op$par$pred]
data[,predi:=opi$par$pred]
#compare to gam
#fit=gam(data=data, formula = y~s(x), family=poisson())
fit=scam(data=data, formula = y~s(x, k=10, bs="mpd", m=2), family=gaussian())
data[,gam:=fit$fitted.values]
#plot result
ggplot(data)+geom_point(aes(x,y),alpha=0.2)+geom_line(aes(x,f),colour="black")+#ylim(0,40)+
  #geom_line(data=mat, aes(x,value,colour=variable))+
  geom_line(aes(x,gam),colour="blue")+
  geom_line(aes(x,predi),colour="green")+
  geom_line(aes(x,pred),colour="red")#+scale_y_log10()






