library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

#matg[,value.ref:=csnorm:::gfl_get_value(valuehat, weight, trails, .1, 1, 0)]
#save(trails,matg,file="tmp_matg.RData")
load("tmp_matg.RData")
ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value.ref))+geom_raster(aes(bin2,bin1,fill=valuehat))+scale_fill_gradient2()

matg[,value:=csnorm:::gfl_get_value(valuehat, weight, trails, .1, 1, 0, nthreads=10)]
matg[,all(value==value.ref)]
matg[,mean(abs(value-value.ref))]
ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+geom_raster(aes(bin2,bin1,fill=value.ref))+scale_fill_gradient2()

