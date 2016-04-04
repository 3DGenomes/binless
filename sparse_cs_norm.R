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

biases=fread("data/rao_HICall_chr19_35000000-36000000_biases.dat")
setkey(biases,id)
counts=fread("data/rao_HICall_chr19_35000000-36000000_counts.dat")

biases=fread("data/caulo_3000000-4000000_biases.dat")
setkey(biases,id)
counts=fread("data/caulo_3000000-4000000_counts.dat")
both=fread("data/caulo_3000000-4000000_both.dat")



#fit it with stan and gam
sm = stan_model(file = "sparse_cs_norm_fit.stan")

system.time(op <- optimizing(sm, data = list(Krow=150, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                             danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                             Kdiag=10, N=counts[,.N],
                                             counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                             cidx=t(data.matrix(counts[,.(id1,id2)]))),
                             as_vector=F, hessian=F, iter=2000, verbose=T, init_alpha=1))
#predict fits
smp = stan_model(file="sparse_cs_norm_predict.stan")
system.time(pred <- optimizing(smp, data = list(Krow=150, Kdiag=10, cutsites=biases[,c(min(pos),max(pos))],
                                                N=10000, S=10000, intercept=op$par$intercept,
                                                eRJ=op$par$eRJ, eDE=op$par$eDE, beta_nu=op$par$beta_nu,
                                                beta_delta=op$par$beta_delta, beta_diag=op$par$beta_diag),
                             as_vector=F, hessian=F, iter=1, verbose=T, init_alpha=1))
biases.pred=data.table(pos=pred$par$genome, nu=exp(pred$par$log_nu), delta=exp(pred$par$log_delta),
                       mean_DL=pred$par$mean_DL, mean_DR=pred$par$mean_DR, mean_RJ=pred$par$mean_RJ)
counts.pred=data.table(distance=pred$par$distnce, decay=exp(pred$par$log_decay))

#compare with previous 6-cutter model
fit=scam(N ~ s(log(distance), bs="mpd", m=2, k=10) + category:(log(dangling.L.1+1) + log(dangling.R.1+1) + log(rejoined.1+1) +
                                         log(dangling.L.2+1) + log(dangling.R.2+1) + log(rejoined.2+1)),
        data=both, family=negbin(1.9))
counts[,decay.gam:=exp(predict(fit, both[category=="contact.up"], type = "terms", terms="s(log(distance))"))]
counts[,mean_cup.gam:=predict(fit, both[category=="contact.up"], type = "response")]

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



ggplot(counts)+scale_y_log10() + geom_point(aes(pos2-pos1,contact.up/mean_cup1*decay1), alpha=0.01) +
  geom_line(aes(pos2-pos1,decay1), colour="red") + geom_line(aes(pos2-pos1,decay.gam),colour="green")
ggplot(counts)+scale_y_log10() + geom_point(aes(pos2-pos1,contact.up/mean_cup.gam*decay.gam), alpha=0.01) +
  geom_line(aes(pos2-pos1,decay1), colour="red") + geom_line(aes(pos2-pos1,decay.gam),colour="green")

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
ggplot(biases) + scale_y_log10() +
  geom_point(aes(pos,dangling.L),colour="red") +
  geom_point(aes(pos,dangling.R),colour="green") +
  geom_point(aes(pos,rejoined),colour="blue") +
  geom_line(data=biases.pred[mean_DL<biases[,max(dangling.L)]&mean_DL>.1], aes(pos,mean_DL),colour="red")+
  geom_line(data=biases.pred[mean_DR<biases[,max(dangling.R)]&mean_DR>.1], aes(pos,mean_DR),colour="green")+
  geom_line(data=biases.pred[mean_RJ<biases[,max(rejoined)]&mean_RJ>.1], aes(pos,mean_RJ),colour="blue")#+xlim(35400000,35415000)


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



