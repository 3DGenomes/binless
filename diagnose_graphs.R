library(csnorm)
library(data.table)
library(ggplot2)
library(foreach)

#setwd("/home/yannick/simulations/cs_norm")
setwd("/Users/yannick/Documents/simulations/cs_norm")

nbins=5

g.R = csnorm:::compute_2d_connectivity_graph(nbins)
plot(g.R)
data=CJ(bin1=1:nbins,bin2=1:nbins)[bin2>=bin1]
data[,value:=(bin1+bin2)%/%2]
ggplot(data)+geom_raster(aes(bin1,bin2,fill=value))

csnorm:::boost_build_patch_graph_components(nbins,data[order(bin2)],1e-3)

a=csnorm:::boost_build_patch_graph_components(nbins,data,1e-3)
b=csnorm:::build_patch_graph(data,csnorm:::gfl_compute_trails(nbins),1e-3)
all(a$membership==b$components$membership-1)

load("tmp_csig.RData")
bic_r = csnorm:::gfl_BIC(csig, .1, 10, -0.1)
ggplot(bic_r$mat)+geom_raster(aes(bin1,bin2,fill=valuehat))+geom_raster(aes(bin2,bin1,fill=value))+scale_fill_gradient2()

bic_c = csnorm:::get_bic(csig, .1, 10, -0.1)

bic_c$dof
bic_r$dof

bic_c$BIC
bic_r$BIC

all.equal(as.data.table(bic_c$mat),bic_r$mat)
a=merge(as.data.table(bic_c$mat)[,.(bin1,bin2,value,patchno=patchno+1)],
        bic_r$mat[,.(bin1,bin2,value,patchno)],by=c("bin1","bin2"),suffixes=c(".c",".r"))




b=foreach (lambda2=5:15,.combine=rbind) %do% {
  bic_r = csnorm:::gfl_BIC(csig, .1, lambda2, -0.1)
  bic_c = csnorm:::get_bic(csig, .1, lambda2, -0.1)
  data.table(lambda2=lambda2,bic.c=bic_c$BIC, bic.r=bic_r$BIC, dof.c=bic_c$dof, dof.r=bic_r$dof)
}
ggplot(b)+geom_line(aes(lambda2,dof.r,colour="r"))+geom_line(aes(lambda2,dof.c,colour="c"))
ggplot(b)+geom_line(aes(lambda2,bic.r,colour="r"))+geom_line(aes(lambda2,bic.c,colour="c"))
