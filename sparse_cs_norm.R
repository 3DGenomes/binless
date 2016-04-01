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
sm = stan_model(file = "sparse_cs_norm_pred.stan")

system.time(op <- optimizing(sm, data = list(Krow=100, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                             danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                             Kdiag=10, N=counts[,.N],
                                             counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                             cidx=t(data.matrix(counts[,.(id1,id2)]))),
                             as_vector=F, hessian=F, iter=1000, verbose=T, init_alpha=1))
#
biases[,mean_RJ1:=op$par$mean_RJ]
biases[,mean_DL1:=op$par$mean_DL]
biases[,mean_DR1:=op$par$mean_DR]
counts[,mean_cclose1:=op$par$mean_cclose]
counts[,mean_cfar1:=op$par$mean_cfar]
counts[,mean_cdown1:=op$par$mean_cdown]
counts[,mean_cup1:=op$par$mean_cup]
counts[,decay1:=exp(op$par$log_decay)]
#
biases[,mean_RJ2:=op$par$mean_RJ]
biases[,mean_DL2:=op$par$mean_DL]
biases[,mean_DR2:=op$par$mean_DR]
counts[,mean_cclose2:=op$par$mean_cclose]
counts[,mean_cfar2:=op$par$mean_cfar]
counts[,mean_cdown2:=op$par$mean_cdown]
counts[,mean_cup2:=op$par$mean_cup]
counts[,decay2:=exp(op$par$log_decay)]
#
biases[,mean_RJ3:=op$par$mean_RJ]
biases[,mean_DL3:=op$par$mean_DL]
biases[,mean_DR3:=op$par$mean_DR]
counts[,mean_cclose3:=op$par$mean_cclose]
counts[,mean_cfar3:=op$par$mean_cfar]
counts[,mean_cdown3:=op$par$mean_cdown]
counts[,mean_cup3:=op$par$mean_cup]
counts[,decay3:=exp(op$par$log_decay)]



ggplot(counts)+scale_y_log10() + geom_point(aes(pos2-pos1,contact.up-mean_cup1)) +
  geom_line(aes(pos2-pos1,decay1), colour="red")
ggplot(counts)+scale_y_log10() + geom_point(aes(pos2-pos1,contact.up-mean_cup2)) +
  geom_line(aes(pos2-pos1,decay2), colour="red")
ggplot(counts)+scale_y_log10() + geom_point(aes(pos2-pos1,contact.up-mean_cup3)) +
  geom_line(aes(pos2-pos1,decay3), colour="red")

ggplot(biases)+scale_y_log10() + geom_point(aes(pos,dangling.L)) +
  geom_line(aes(pos,mean_DL1))+
  geom_line(aes(pos,mean_DL2))+
  geom_line(aes(pos,mean_DL3))#+xlim(35400000,35405000)
ggplot(biases)+scale_y_log10() + geom_point(aes(pos,dangling.R)) +
  geom_line(aes(pos,mean_DR1))+
  geom_line(aes(pos,mean_DR2))+
  geom_line(aes(pos,mean_DR3))
ggplot(biases)+scale_y_log10() + geom_point(aes(pos,rejoined)) +
  geom_line(aes(pos,mean_RJ1))+
  geom_line(aes(pos,mean_RJ2))+
  geom_line(aes(pos,mean_RJ3))
ggplot(biases)+scale_y_log10() +
  geom_point(aes(pos,dangling.L),colour="red") +
  geom_point(aes(pos,dangling.R),colour="green") +
  geom_point(aes(pos,rejoined),colour="blue") +
  geom_line(aes(pos,mean_DL1),colour="red")+
  geom_line(aes(pos,mean_DR1),colour="green")+
  geom_line(aes(pos,mean_RJ1),colour="blue")+xlim(35400000,35405000)


#data[,delta:=exp(op$par$log_delta)]
#data[,fij:=op$par$decay]
#plot result
ggplot()+#scale_y_log10()+
  geom_line(data=melt(biases,id.vars=c("id","pos"))[variable%in%c("nu1","nu2","nu3")],
            aes(pos,value,colour=variable))#+
  geom_point(data=melt(biases,id.vars=c("id","pos"))[!(variable%in%c("nu1","nu2","nu3"))],
            aes(pos,value,colour=variable))

ggplot()+#scale_y_log10()+
    geom_line(data=melt(biases,id.vars=c("id","pos"))[variable%in%c("delta1","delta2","delta3")],
              aes(pos,value,colour=variable))




            aes(pos2-pos1,value,colour=variable))+
  geom_point(data=melt(counts,id.vars=c("id1","id2","pos1","pos2"))[!(variable%in%c("fij1","fij2","fij3"))],
            aes(pos2-pos1,value,colour=variable))

  
ggplot(melt(data,id.vars=c("id","pos")),aes(pos,value,colour=variable))+geom_point()+geom_line()+scale_y_log10()
ggplot(melt(data,id.vars=c("id","pos")),aes(pos,value,colour=variable))+geom_point()+geom_line()+scale_y_log10()+xlim(35450000, 35475000)
ggplot(data) +
  geom_line(aes(pos,nu),colour="green")+
  geom_line(aes(pos,delta),colour="blue")



