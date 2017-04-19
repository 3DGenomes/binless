library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(microbenchmark)

setwd("/home/yannick/simulations/cs_norm")

#matg[,value.ref:=csnorm:::gfl_get_value(valuehat, weight, trails, .1, 1, 0)]
#save(trails,matg,file="tmp_matg.RData")
load("tmp_matg.RData")
ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value.ref))+geom_raster(aes(bin2,bin1,fill=valuehat))+scale_fill_gradient2()

matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, .1, 1, 0, nthreads=30)]
matg[,all(value==value.ref)]
matg[,mean(abs(value-value.ref))]
ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+geom_raster(aes(bin2,bin1,fill=value.ref))+scale_fill_gradient2()

bench=function(x){matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, .1, 1, 0, nthreads=x)]}
jobs = lapply(1:17, function(s) local({s = s; bquote(bench(.(s)))}) )
a=microbenchmark(list=jobs,times=100) 
plot(a)
a.dt=data.table(call=a$expr,time=a$time)[,.(call,time,ncpus=unclass(call))]
a.dt[,speedup:=.SD[ncpus==1,mean(time)]/time]
ggplot(a.dt)+geom_jitter(aes(ncpus,speedup))
ggplot(a.dt)+geom_boxplot(aes(ncpus,speedup,group=ncpus))
ggplot(a.dt)+geom_boxplot(aes(ncpus,speedup/ncpus,group=ncpus))


csnorm:::gfl_BIC(matg, trails, .1, 1, 0, tol.value=1e-3, ncores=30)

lambda2 = csnorm:::optimize_lambda2(matg, trails, tol=1e-3, ncores=30)
profvis({vals = csnorm:::optimize_lambda1_eCprime(matg, trails, tol=1e-3, lambda2=1, positive=T, ncores=1)})
vals = csnorm:::optimize_lambda1_only(matg, trails, tol=1e-3, lambda2=1, positive=T, ncores=30)
