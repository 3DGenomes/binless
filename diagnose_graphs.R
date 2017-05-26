library(csnorm)
library(data.table)
library(ggplot2)

setwd("/home/yannick/simulations/cs_norm")

nbins=5

g.R = csnorm:::compute_2d_connectivity_graph(nbins)
plot(g.R)
data=CJ(bin1=1:nbins,bin2=1:nbins)[bin2>=bin1]
data[,value:=(bin1+bin2)%/%2]
ggplot(data)+geom_raster(aes(bin1,bin2,fill=value))

csnorm:::boost_build_patch_graph(nbins,data[order(bin2)],1e-3)

a=csnorm:::boost_build_patch_graph(nbins,data,1e-3)
b=csnorm:::build_patch_graph(data,csnorm:::gfl_compute_trails(nbins),1e-3)

