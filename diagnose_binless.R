library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)


setwd("/home/yannick/simulations/cs_norm")

load("data/cohesin_MboI_single_csnorm_optimized.RData")

dispersion=cs@par$alpha
trails=cs@settings$trails
tol.val=cs@settings$tol.val
lambda2=10

cts = csnorm:::csnorm_gauss_signal_muhat_mean(cs, cs@zeros, cs@settings$sbins)[
  name==name[1],.(name,bin1,bin2,count,lmu.nosig,phi,weight)]
setkeyv(cts,c("name","bin1","bin2"))
trails$nrow = cts[,max(unclass(bin2))]

### performance iteration: R
niter=2
values=data.table()
cts[,phi:=0]
for (i in 1:niter) {
  #compute mat
  cts[,mu:=exp(lmu.nosig+phi)]
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/dispersion))]
  
  mat = cts[,.(phihat=weighted.mean(z+phi, weight/var),
               phihat.var=2/sum(weight/var),
               ncounts=sum(weight)),keyby=c("name","bin1","bin2")]
  mat = mat[cs@par$signal[,.(name,bin1,bin2)],,on=c("name","bin1","bin2")] #to add empty rows/cols
  mat[is.na(phihat),c("phihat","phihat.var","ncounts"):=list(1,Inf,0)] #bins with no detectable counts
  mat[,c("valuehat","weight","diag.idx"):=list(phihat,1/phihat.var,unclass(bin2)-unclass(bin1))]
  stopifnot(mat[weight==0,.N==0])
  #mat[diag.idx<=1,weight:=0]
  setkey(mat,name,bin1,bin2)
  
  #GFL
  alpha=0.2
  inflate=2
  maxsteps=10000
  value.old=value
  value = csnorm:::weighted_graphfl(mat[,valuehat], mat[,weight], trails$ntrails, trails$trails,
                                    trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.val/2)
  cat("epsilon=",mean(abs(value-value.old)))
  
  #report phi
  mat[,phi:=value]
  values=rbind(values,mat[,.(step=i,bin1,bin2,phi)])
  cts[,phi:=NULL]
  cts=merge(cts,mat[,.(name,bin1,bin2,phi)],keyby=c("name","bin1","bin2"))
}
value.r=value

ggplot(values)+geom_raster(aes(bin1,bin2,fill=phi))+facet_wrap(~step)+scale_fill_gradient2()



### performance iteration: C++

niter=10
cts[,phi:=0]
alpha=0.2
inflate=2
maxsteps=10000


value.c=csnorm:::wgfl_perf(cts, cs@par$alpha, niter, trails$nrow, trails$ntrails, trails$trails,
                           trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.val/2)

all.equal(value.c,value.r)



