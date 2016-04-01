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

biases=fread("caulo_toy_biases.dat")
counts=fread("caulo_toy_counts.dat")
data=unique(rbind(counts[,.(id=id1,pos=pos1, rejoined=rejoined.1)],counts[,.(id=id2,pos=pos2)]))[biases,,on="id"]


data=fread("cs_norm_caulo.dat", stringsAsFactors = T)
counts=fread("cs_norm_rao_HiC035_chr22.dat", stringsAsFactors = T)
data=unique(rbind(counts[,.(id=id1,pos=pos1)],counts[,.(id=id2,pos=pos2)]))[biases,,on="id"]
biases=fread("rao_HICall_chr19_35400000-35500000_biases.dat")
setkey(biases,id)

data=biases
ggplot(melt(data,id.vars=c("id","pos")),aes(pos,value,colour=variable))+geom_point()+geom_line()+scale_y_log10()

counts=fread("rao_HICall_chr19_35400000-35500000_counts.dat")

#fit it with stan and gam
sm = stan_model(file = "sparse_cs_norm.stan")

system.time(op <- optimizing(sm, data = list(Krow=5, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                             danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                             Kdiag=5, N=counts[,.N],
                                             counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                             cidx=t(data.matrix(counts[,.(id1,id2)]))),
                             as_vector=F, hessian=F, iter=1000, verbose=T))
biases[,nu1:=exp(op$par$log_nu)]
biases[,delta1:=exp(op$par$log_delta)]
counts[,fij1:=exp(op$par$log_decay)]

biases[,nu2:=exp(op$par$log_nu)]
biases[,delta3:=exp(op$par$log_delta)]
counts[,fij2:=exp(op$par$log_decay)]

biases[,nu3:=exp(op$par$log_nu)]
biases[,delta3:=exp(op$par$log_delta)]
counts[,fij3:=exp(op$par$log_decay)]

#data[,delta:=exp(op$par$log_delta)]
#data[,fij:=op$par$decay]
#plot result
ggplot()+scale_y_log10()+
  geom_line(data=melt(data,id.vars=c("id","pos"))[variable%in%c("nu1","nu2","nu3")],
            aes(pos,value,colour=variable))+
  geom_point(data=melt(data,id.vars=c("id","pos"))[!(variable%in%c("nu1","nu2","nu3"))],
            aes(pos,value,colour=variable))


ggplot(melt(data,id.vars=c("id","pos")),aes(pos,value,colour=variable))+geom_point()+geom_line()+scale_y_log10()
ggplot(melt(data,id.vars=c("id","pos")),aes(pos,value,colour=variable))+geom_point()+geom_line()+scale_y_log10()+xlim(35450000, 35475000)
ggplot(data) +
  geom_line(aes(pos,nu),colour="green")+
  geom_line(aes(pos,delta),colour="blue")



