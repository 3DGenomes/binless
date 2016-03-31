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

setnames(biases, "re.pos","pos")
data=biases
ggplot(melt(data,id.vars=c("id","pos")),aes(pos,value,colour=variable))+geom_point()+geom_line()+scale_y_log10()

counts=fread("rao_HICall_chr19_35400000-35500000_counts.dat")
setnames(counts, "re.pos1","pos1")
setnames(counts, "re.pos2","pos2")
counts=dcast(counts[,.(category,pos1,pos2,N)], pos1+pos2~category, fill=0, value.var="N")

#fit it with stan and gam
sm = stan_model(file = "sparse_cs_norm.stan")

system.time(op <- optimizing(sm, data = list(Krow=5, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                             danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                             Kdiag=5, N=counts[,.N], ),
                             
                             as_vector=F, hessian=F, iter=1000, verbose=F))
data[,nu1:=exp(op$par$log_nu)]
system.time(op <- optimizing(sm, data = list(Krow=5, S=data[,.N], cutsites=data[,pos], rejoined=data[,rejoined],
                                             danglingL=data[,dangling.L], danglingR=data[,dangling.R]),
                             as_vector=F, hessian=F, iter=1000, verbose=F))
data[,nu2:=exp(op$par$log_nu)]
system.time(op <- optimizing(sm, data = list(Krow=5, S=data[,.N], cutsites=data[,pos], rejoined=data[,rejoined],
                                             danglingL=data[,dangling.L], danglingR=data[,dangling.R]),
                             as_vector=F, hessian=F, iter=1000, verbose=F))
data[,nu3:=exp(op$par$log_nu)]


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



