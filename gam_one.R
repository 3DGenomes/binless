library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)

setwd("/home/yannick/simulations/cs_norm")

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
#data=generate_data(fun=function(x){3+2*sin(5*x)+3*x}, sd=100, xmin=0,xmax=3, npoints=1000)
#data[,pos:=x]
#data[,x:=x*1e5/3.+35400000]
#data=generate_data(fun=function(x){1+5*exp(-(x-1.5)**2*2)}, sd=5, xmin=0,xmax=3, npoints=1000)
#data=data[y>0]
#data=generate_data(fun=function(x){5*sin(3*x)+5.1}, sd=1, xmin=0,xmax=3, npoints=500)
data=generate_data(fun=function(x){300+20*sin(x)+.1*x}, sd=10, xmin=0,xmax=10, npoints=100)
ggplot(data)+geom_point(aes(x,y))+geom_line(aes(x,f),colour="red") 


#fit it with stan
sm = stan_model(file = "gam_one_sparse.stan")
op = optimizing(sm, data = list(N=data[,.N], K=5, y=data[,y], x=data[,x]),
                as_vector=F, hessian=F, iter=100000, verbose=T)
data[,pred:=op$par$pred]
fit = gam(data=data, formula = y ~ s(x,bs="ps", m=2, k=5,by=one))
data[,gam:=fit$fitted.values]
ggplot(data)+geom_line(aes(x,gam,colour="gam"))+geom_line(aes(x,pred,colour="pred"))
#
mat.pred=stan_matrix_to_datatable(op$par$designmat, data[,x])[,.(x,variable,value,method="pred")]
mat.gam=stan_matrix_to_datatable(predict(fit,type="lpmatrix"), data[,x])[,.(x,variable,value,method="gam")]
mat=rbindlist(list(pred=mat.pred,gam=mat.gam), use.names = T)
ggplot(mat)+geom_line(aes(x,value,colour=variable))+facet_wrap(~method)









load("biases.RData")
data=rbind(biases[cat%in%c("dangling L","contact L"),
                  .(cat=factor(cat),x=pos-min(pos),y=log_iota+z,std,weight=1/std^2,iota.coef=1,rho.coef=0)],
           biases[cat%in%c("dangling R","contact R"),
                  .(cat=factor(cat),x=pos-min(pos),y=log_rho+z,std,weight=1/std^2, iota.coef=0,rho.coef=1)],
           biases[cat=="rejoined",
                  .(cat=factor(cat),x=pos-min(pos),y=(log_iota+log_rho)/2+z,std,weight=1/std^2, iota.coef=1/2, rho.coef=1/2)])
data[,is.rejoined:=ifelse(cat=="rejoined",1,0)]
data[,is.dangling:=ifelse(cat %in% c("dangling L","dangling R"),1,0)]
data[cat == "rejoined", y:=y+1]
data[cat %in% c("dangling L","dangling R"), y:=y+5]
data[cat %in% c("contact L","contact R"), y:=y+10]
data=data[x<1e5]
data=rbind(cbind(data,dset="a"),cbind(data,dset="b"))
data[dset=="b",y:=y+10]
data[,dset:=as.integer(factor(dset))-1]
ggplot(data)+geom_pointrange(aes(x,y,ymin=y-std,ymax=y+std,colour=cat),alpha=0.1)+facet_wrap(~dset)
ggplot(data)+geom_pointrange(aes(x,y,ymin=y-std,ymax=y+std,colour=cat))+xlim(0,3e4)




#toy fit
fit=gam(data=data[,.SD[1:5],by=c("cat","dset")], formula = y ~ dset+is.dangling+is.rejoined-1+s(x,bs="ps", m=2, k=4, by=iota.coef)+s(x,bs="ps", m=2, k=4, by=rho.coef),
        family=gaussian(), scale=1, fit=F)
fit=gam(data=data[,.SD[1:5],by=cat], formula = y ~ is.rejoined-1 + s(x, by=iota.coef,bs="ps", m=2, k=4) + s(x, by=rho.coef,bs="ps", m=2, k=4),
        family=gaussian(), scale=1, fit=F)
#X1=predict(fit, type="lpmatrix")
colSums(fit$X)
dim(fit$X)
rankMatrix(fit$X)
a=data.table(fit$X)
a[,idx:=.I]
a=melt(a, id.vars="idx")
ggplot(a)+geom_line(aes(idx,value,colour=variable))+xlim(0,5)

#fit genomic biases
fit=gam(data=data, formula = y ~ dset+is.dangling+is.rejoined-1 + s(x, by=iota.coef,bs="ps", m=2, k=50) + s(x, by=rho.coef,bs="ps", m=2, k=50),
        family=gaussian(), weight=weight, scale=1)
data[,gam:=fit$fitted.values]
dcast(data[,.(cat,x,gam,y,dset)], dset+x ~ cat, value.var=c("gam","y"))
ggplot(data)+geom_pointrange(aes(x,y,ymin=y-std,ymax=y+std,colour=cat),alpha=0.1)+
  geom_line(aes(x,gam,colour=cat))+facet_wrap(~dset)


#plot result
ggplot(data[cat=="contact L"])+#geom_point(aes(x,y),alpha=0.2)+#geom_line(aes(x,f),colour="black")+#ylim(0,40)+
  geom_pointrange(aes(x,y,ymin=y-std,ymax=y+std),alpha=0.1)+
  geom_line(aes(x,gam),colour="blue") + xlim(0,1e5) #+
#geom_line(aes(x,pred),colour="red") #+scale_y_log10()
ggplot(data[cat=="dangling L"])+#geom_point(aes(x,y),alpha=0.2)+#geom_line(aes(x,f),colour="black")+#ylim(0,40)+
  geom_pointrange(aes(x,y,ymin=y-std,ymax=y+std),alpha=0.1)+
  geom_line(aes(x,gam),colour="blue")+xlim(0,1e5)

