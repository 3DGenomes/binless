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

csnorm:::boost_get_number_of_patches(nbins,data,1e-3)
