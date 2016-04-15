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

stan_matrix_to_datatable = function(opt, x) {
  vals=data.table(opt)
  vals[,x:=x]
  melt(data.table(vals), id.vars="x")
}

optimize_exposures = function(model, biases, counts, iter=1000000, verbose=T) {
  optimizing(model, data = list(S=biases[,.N], rejoined=biases[,rejoined],
                                danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                N=counts[,.N],
                                counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)]))),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

optimize_nu = function(model, biases, op0, Krow=1000, lambda_nu=1, iter=1000000, verbose=T) {
  optimizing(model, data = list(Krow=Krow, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                eRJ=op0$par$eRJ, eDE=op0$par$eDE, lambda_nu=lambda_nu),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

optimize_delta = function(model, biases, op0, op1, Krow=1000, lambda_delta=1, iter=1000000, verbose=T) {
  optimizing(model, data = list(Krow=Krow, S=biases[,.N], cutsites=biases[,pos],
                               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                               eDE=op0$par$eDE, beta_nu=op1$par$beta_nu, lambda_delta=lambda_delta),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

optimize_decay = function(model, biases, counts, op0, op1, op2, Kdiag=10, lambda_diag=1, iter=1000000, verbose=T) {
  optimizing(model, data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], N=counts[,.N],
                                 counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                 cidx=t(data.matrix(counts[,.(id1,id2)])),
                                 eC=op0$par$eC, eRJ=op0$par$eRJ, eDE=op0$par$eDE,
                                 log_nu=op1$par$log_nu, log_delta=op2$par$log_delta, lambda_diag=lambda_diag),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

optimize_all = function(model, biases, counts, op0, op1, op2, op3, Krow=1000, Kdiag=10,
                        lambda_nu=1, lambda_delta=1, lambda_diag=1, iter=1000000, verbose=T) {
  optimizing(model, data = list( Krow=Krow, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                 danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                 Kdiag=Kdiag, N=counts[,.N],
                                 counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                 cidx=t(data.matrix(counts[,.(id1,id2)])),
                                 lambda_nu=lambda_nu, lambda_delta=lambda_delta, lambda_diag=lambda_diag),
             init = list(eC=op0$par$eC, eRJ=op0$par$eRJ, eDE=op0$par$eDE,
                         beta_nu=op1$par$beta_nu, beta_delta=op2$par$beta_delta, beta_diag=op3$par$beta_diag,
                         alpha=op3$par$alpha),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}

optimize_all_noinit = function(model, biases, counts, Krow=1000, Kdiag=10, 
                               lambda_nu=1, lambda_delta=1, lambda_diag=1, iter=1000000, verbose=T) {
  optimizing(model, data = list( Krow=Krow, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                 danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                 Kdiag=Kdiag, N=counts[,.N],
                                 counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                 cidx=t(data.matrix(counts[,.(id1,id2)])),
                                 lambda_nu=lambda_nu, lambda_delta=lambda_delta, lambda_diag=lambda_diag),
             as_vector=F, hessian=F, iter=iter, verbose=verbose)
}






biases=fread("data/rao_HICall_chr20_all_biases.dat")
setkey(biases,id)
counts=fread("data/rao_HICall_chr20_all_counts.dat")
both=fread("data/rao_HICall_chr20_all_both.dat")


biases=fread("data/rao_HIC035_chr20_all_biases.dat")
setkey(biases,id)
counts=fread("data/rao_HIC035_chr20_all_counts.dat")
both=fread("data/rao_HIC035_chr20_all_both.dat")


biases=fread("data/rao_HICall_chr20_35000000-36000000_biases.dat")
setkey(biases,id)
counts=fread("data/rao_HICall_chr20_35000000-36000000_counts.dat")
both=fread("data/rao_HICall_chr20_35000000-36000000_both.dat")

biases=fread("data/rao_HIC035_chr20_35000000-36000000_biases.dat")
setkey(biases,id)
counts=fread("data/rao_HIC035_chr20_35000000-36000000_counts.dat")
both=fread("data/rao_HIC035_chr20_35000000-36000000_both.dat")



biases=fread("data/caulo_3000000-4000000_biases.dat")
setkey(biases,id)
counts=fread("data/caulo_3000000-4000000_counts.dat")
both=fread("data/caulo_3000000-4000000_both.dat")

biases=fread("data/caulo_all_biases.dat")
setkey(biases,id)
counts=fread("data/caulo_all_counts.dat")
both=fread("data/caulo_all_both.dat")
counts=counts[sample(.N,min(25000,.N))]







#make sure the 6-cutter model fits it nicely
sub=both[,.SD[sample(.N,min(25000,.N))],by=category]
fit=gam(N ~ s(log(distance)) + category:(log(dangling.L.1+1) + log(dangling.R.1+1) + log(rejoined.1+1) +
                                           log(dangling.L.2+1) + log(dangling.R.2+1) + log(rejoined.2+1)),
        data=sub, family=nb())
summary(fit)
sub[,fij:=exp(predict.gam(fit, sub, type="terms", terms="s(log(distance))"))]
sub[,mean:=fit$fitted.values]
ggplot(sub[distance>70e4])+
  geom_point(aes(distance,N*fij/mean),alpha=0.01)+
  geom_line(aes(distance,fij))+scale_x_log10()+scale_y_log10()



#### simple optimization
counts=counts[sample(.N,25000)]
smfit = stan_model(file = "sparse_cs_norm_fit.stan")
system.time(op <- optimize_all_noinit(smfit, biases, counts, Krow=10000, Kdiag=10,
                                      lambda_nu=1, lambda_delta=1, lambda_diag=1))
#compare deviances
c(100*(fit$null.deviance-fit$deviance)/fit$null.deviance, op$par$deviance_proportion_explained)
#compare dispersions
c(fit$family$getTheta(T), op$par$alpha)
#compare decay
a=cbind(counts,data.table(decay=exp(op$par$log_decay),
                          cclose=exp(op$par$log_mean_cclose),
                          cfar=exp(op$par$log_mean_cfar),
                          cup=exp(op$par$log_mean_cup),
                          cdown=exp(op$par$log_mean_cdown)))
ggplot(a[pos2-pos1>1e4][sample(.N,min(.N,10000))])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(pos2-pos1, contact.close*decay/cclose), colour="blue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.far*decay/cclose), colour="green", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.up*decay/cclose), colour="darkblue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.down*decay/cclose), colour="darkgreen", alpha=0.01)+
  geom_line(aes(pos2-pos1, decay) ,colour="pink")


smfit = stan_model(file = "sparse_cs_norm_fit.stan")
system.time(op.sep <- optimize_all_noinit(smfit, biases, counts, Krow=10000, Kdiag=10,
                                      lambda_nu=1, lambda_delta=1, lambda_diag=1))
#compare deviances
c(op$par$deviance_proportion_explained, op.sep$par$deviance_proportion_explained,
  op.new$par$deviance_proportion_explained)
#compare dispersions
c(op$par$alpha, op.sep$par$alpha, op.sep$par$aDE, op.sep$par$aRJ)
c(op$par$alpha, op.new$par$aC, op.new$par$aDE, op.new$par$aRJ)
#compare decay
a=cbind(counts,data.table(decay=exp(op$par$log_decay),
                          cclose=exp(op$par$log_mean_cclose)),
        data.table(decay.sep=exp(op.sep$par$log_decay),
                   cclose.sep=exp(op.sep$par$log_mean_cclose)),
        data.table(decay.new=exp(op.new$par$log_decay),
                   cclose.new=exp(op.new$par$log_mean_cclose)))
ggplot(a[pos2-pos1>1e4])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(pos2-pos1, contact.close*decay.new/cclose.new), colour="blue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.far*decay.new/cclose.new), colour="green", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.up*decay.new/cclose.new), colour="darkblue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.down*decay.new/cclose.new), colour="darkgreen", alpha=0.01)+
  geom_line(aes(pos2-pos1, decay) ,colour="pink")+
  geom_line(aes(pos2-pos1, decay.new) ,colour="red")+
  geom_line(aes(pos2-pos1, decay.sep) ,colour="purple")+
  geom_line(aes(distance, fij), colour="blue", data=sub[distance>1e4])






### initialization + optimization
#initial guesses for exposures
smb0 = stan_model(file = "sparse_cs_norm_init_0_exposures.stan")
system.time(op.b0 <- optimize_exposures(smb0, biases, counts))
op2.b0 <- optimize_exposures(smb0, biases, counts)
#
a=data.table(op1=as.numeric(cbind(op.b0$par)[,1]), op2=as.numeric(cbind(op2.b0$par)[,1]),
             name=names(op2.b0$par))
ggplot(a)+geom_bar(aes(name,op2-op1),stat="identity", position="dodge")+ylim(-.1,.1)
ggplot(a)+geom_bar(aes(name,(op2-op1)/max(op2,op1)),stat="identity", position="dodge")+ylim(-.01,.01)


#initial guesses for nu, need to set lambda appropriately
smb1 = stan_model(file = "sparse_cs_norm_init_1_nu.stan")
op.b1 <- optimize_nu(smb1, biases, op.b0, Krow=65000, lambda_nu=1)
op2.b1 <- optimize_nu(smb1, biases, op.b0, Krow=65000, lambda_nu=1)
#predict nu everywhere
pred = stan_model(file = "predict_spline_sparse.stan")
op.pred_nu <- optimizing(pred, data = list(K=65000, N=100000, xrange=biases[,c(min(pos),max(pos))],
                                        intercept=0, beta=op.b1$par$beta_nu),
                      as_vector=F, hessian=F, iter=1, verbose=F)
op2.pred_nu <- optimizing(pred, data = list(K=65000, N=100000, xrange=biases[,c(min(pos),max(pos))],
                                        intercept=0, beta=op2.b1$par$beta_nu),
                      as_vector=F, hessian=F, iter=1, verbose=F)
nu.pred.weighted <- stan_matrix_to_datatable(exp(op.pred_nu$par$weighted), op.pred_nu$par$x)
nu.pred <- data.table(pos=op.pred_nu$par$x, nu1=exp(op.pred_nu$par$log_mean), nu2=exp(op2.pred_nu$par$log_mean))
#
a=data.table(name=c("alpha","deviance_proportion_explained"),
             op1=c(op.b1$par$alpha, op.b1$par$deviance_proportion_explained), 
             op2=c(op2.b1$par$alpha, op2.b1$par$deviance_proportion_explained))
a
ggplot(a)+geom_bar(aes(name,(op2-op1)/max(op2,op1)),stat="identity", position="dodge")#+scale_y_continuous(-.01,.01)
#
a=cbind(biases,data.table(op1=exp(op.b1$par$log_nu), op2=exp(op2.b1$par$log_nu)))
pbegin=35100000 #3166716
pend=35200000 #3191637
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op.b0$par$eDE)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op.b0$par$eDE)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op.b0$par$eRJ)),colour="red")+
  geom_point(aes(pos, op1),colour="blue")+
  geom_point(aes(pos, op2),colour="green")+
  geom_line(data=nu.pred[pos>=pbegin&pos<=pend], aes(pos, nu1), colour="blue")+
  geom_line(data=nu.pred[pos>=pbegin&pos<=pend], aes(pos, nu2), colour="green")#+
  #geom_line(data=nu.pred.weighted[x>=pbegin&x<=pend], aes(x,value,colour=variable))
message("MAD for nu: ", mad(a$op1-a$op2))
#




#initial guesses for delta, need to set lambda appropriately
smb2 = stan_model(file = "sparse_cs_norm_init_2_delta.stan")
op.b2 <- optimize_delta(smb2, biases, op.b0, op.b1, Krow=65000, lambda_delta=1)
op2.b2 <- optimize_delta(smb2, biases, op.b0, op.b1, Krow=65000, lambda_delta=1)
#predict delta everywhere
pred = stan_model(file = "predict_spline_sparse.stan")
op.pred_delta <- optimizing(pred, data = list(K=65000, N=100000, xrange=biases[,c(min(pos),max(pos))],
                                        intercept=0, beta=op.b2$par$beta_delta),
                      as_vector=F, hessian=F, iter=1, verbose=F)
op2.pred_delta <- optimizing(pred, data = list(K=65000, N=100000, xrange=biases[,c(min(pos),max(pos))],
                                         intercept=0, beta=op2.b2$par$beta_delta),
                       as_vector=F, hessian=F, iter=1, verbose=F)
delta.pred.weighted <- stan_matrix_to_datatable(exp(op.pred_delta$par$weighted), op.pred_delta$par$x)
delta.pred <- data.table(pos=op.pred_delta$par$x, delta1=exp(op.pred_delta$par$log_mean), delta2=exp(op2.pred_delta$par$log_mean))
#
data.table(name=c("alpha","deviance_proportion_explained"),
             op1=c(op.b2$par$alpha, op.b2$par$deviance_proportion_explained), 
             op2=c(op2.b2$par$alpha, op2.b2$par$deviance_proportion_explained))
#DE op1 vs op2, centered on RJ
a=cbind(biases,data.table(op1=exp(op.b2$par$log_delta),
                          op2=exp(op2.b2$par$log_delta),
                          log_nu=op.b1$par$log_nu))
pbegin=35100000 #3166716
pend=35200000 #3191637
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/(exp(log_nu+op.b0$par$eDE))),colour="orange")+
  geom_point(aes(pos, dangling.R/(exp(log_nu+op.b0$par$eDE))),colour="pink")+
  geom_point(aes(pos, rejoined/(exp(log_nu+op.b0$par$eRJ))),colour="red")+
  geom_point(aes(pos, op1),colour="blue", shape=0)+
  geom_point(aes(pos, op2),colour="green", shape=0)+
  geom_point(aes(pos, 1/op1),colour="blue", shape=0)+
  geom_point(aes(pos, 1/op2),colour="green", shape=0)+
  geom_line(data=delta.pred[pos>=pbegin&pos<=pend], aes(pos, delta1), colour="blue")+
  geom_line(data=delta.pred[pos>=pbegin&pos<=pend], aes(pos, delta2), colour="green")+
  geom_line(data=delta.pred[pos>=pbegin&pos<=pend], aes(pos, 1/delta1), colour="blue", linetype=2)+
  geom_line(data=delta.pred[pos>=pbegin&pos<=pend], aes(pos, 1/delta2), colour="green", linetype=2)+
  geom_line(data=delta.pred.weighted[x>=pbegin&x<=pend], aes(x,value,colour=variable))
#all 3 types
a=cbind(biases,data.table(RJ=exp(op.b0$par$eRJ+op.b1$par$log_nu),
                          DL=exp(op.b0$par$eDE+op.b1$par$log_nu+op.b2$par$log_delta),
                          DR=exp(op.b0$par$eDE+op.b1$par$log_nu-op.b2$par$log_delta)))
a.pred=data.table(pos=nu.pred$pos, RJ=exp(op.b0$par$eRJ)*nu.pred$nu1,
                  DL=exp(op.b0$par$eDE)*nu.pred$nu1*delta.pred$delta1,
                  DR=exp(op.b0$par$eDE)*nu.pred$nu1/delta.pred$delta1)
ggplot(a[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, rejoined), colour="red")+ geom_point(aes(pos, RJ), shape=0, colour="red")+ geom_line(data=a.pred[pos>=pbegin&pos<=pend], aes(pos, RJ), colour="red")+
  geom_point(aes(pos, dangling.L), colour="orange")+ geom_point(aes(pos, DL), shape=0, colour="orange")+ geom_line(data=a.pred[pos>=pbegin&pos<=pend], aes(pos, DL), colour="orange")+
  geom_point(aes(pos, dangling.R), colour="pink")+ geom_point(aes(pos, DR), shape=0, colour="pink")+ geom_line(data=a.pred[pos>=pbegin&pos<=pend], aes(pos, DR), colour="pink")


#initial guesses for diagonal decay, need to set lambda appropriately
smb3 = stan_model(file = "sparse_cs_norm_init_3_decay.stan")
op.b3 <- optimize_decay(smb3, biases, counts, op.b0, op.b1, op.b2, Kdiag=20, lambda_diag=.1)
op2.b3 <- optimize_decay(smb3, biases, counts, op.b0, op.b1, op.b2, Kdiag=20, lambda_diag=.1)
#
data.table(name=c("alpha","deviance_proportion_explained"),
             op1=c(op.b3$par$alpha, op.b3$par$deviance_proportion_explained), 
             op2=c(op2.b3$par$alpha, op2.b3$par$deviance_proportion_explained))
#DE op1 vs op2, centered on RJ
a=cbind(counts,data.table(decay1=exp(op.b3$par$log_decay),
                          decay2=exp(op2.b3$par$log_decay),
                          cclose1=exp(op.b3$par$log_mean_cclose),
                          cclose2=exp(op2.b3$par$log_mean_cclose),
                          cfar1=exp(op.b3$par$log_mean_cfar),
                          cfar2=exp(op2.b3$par$log_mean_cfar),
                          cup1=exp(op.b3$par$log_mean_cup),
                          cup2=exp(op2.b3$par$log_mean_cup),
                          cdown1=exp(op.b3$par$log_mean_cdown),
                          cdown2=exp(op2.b3$par$log_mean_cdown)))
ggplot(a[pos2-pos1>1e4][sample(.N,min(.N,10000))])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(pos2-pos1, contact.close*decay1/cclose1), colour="blue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.far*decay1/cclose1), colour="green", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.up*decay1/cclose1), colour="darkblue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.down*decay1/cclose1), colour="darkgreen", alpha=0.01)+
  geom_line(aes(pos2-pos1, decay1) ,colour="pink")+
  geom_line(aes(pos2-pos1, decay2) ,colour="orange")

ggplot(counts[pos2-pos1>1e4])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(pos2-pos1, contact.close), colour="blue", alpha=0.01)

counts[,dbin.100:=cut2(log(pos2-pos1), g = 100, levels.mean=T)]
a.100=counts[pos2-pos1>1e4][,sum(contact.close),by=as.numeric(as.character(dbin.100))]
setnames(a.100, c("dist","count"))
counts[,dbin.1000:=cut2(log(pos2-pos1), g = 1000, levels.mean=T)]
a.1000=counts[pos2-pos1>1e4][,sum(contact.close),by=as.numeric(as.character(dbin.1000))]
setnames(a.1000, c("dist","count"))
ggplot(a.1000)+scale_y_log10()+
  geom_point(data=a.1000, aes(dist,count), colour="blue")+
  geom_point(data=a.100, aes(dist,count), colour="red")+
  geom_point(data=counts[sample(.N,10000)][pos2-pos1>1e4], aes(log(pos2-pos1), contact.close), colour="green", alpha=0.1)


fit=gam(count ~ s(dist), data=a.1000, family=nb())
summary(fit)
a.1000[,fij:=fit$fitted.values]
ggplot(a.1000)+
  #geom_point(aes(distance,count))+
  geom_line(aes(distance,fij))+scale_x_log10()+scale_y_log10()
lm(data=a.1000[distance<7e5], log(fij)~log(distance))




###final optimization
smfit = stan_model(file = "sparse_cs_norm_fit.stan")
system.time(op.all <- optimize_all(smfit, biases, counts, op.b0, op.b1, op.b2, op.b3, Krow=1000, Kdiag=10,
                      lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
system.time(op2.all <- optimize_all(smfit, biases, counts, Krow=1000, Kdiag=10,
                       lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
#compare deviances
data.table(name=c("alpha","deviance_proportion_explained"),
           op1=c(op.all$par$alpha, op.all$par$deviance_proportion_explained), 
           op2=c(op2.all$par$alpha, op2.all$par$deviance_proportion_explained))
#predict nu everywhere
pred = stan_model(file = "predict_spline_sparse.stan")
op.all.pred_nu <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                           intercept=0, beta=op.all$par$beta_nu),
                         as_vector=F, hessian=F, iter=1, verbose=F)
op2.all.pred_nu <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                            intercept=0, beta=op2.all$par$beta_nu),
                          as_vector=F, hessian=F, iter=1, verbose=F)
nu.all.pred.weighted <- stan_matrix_to_datatable(exp(op.all.pred_nu$par$weighted), op.all.pred_nu$par$x)
nu.all.pred <- data.table(pos=op.all.pred_nu$par$x, nu1=exp(op.all.pred_nu$par$log_mean), nu2=exp(op2.all.pred_nu$par$log_mean))
#predict delta everywhere
op.all.pred_delta <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                              intercept=0, beta=op.all$par$beta_delta),
                            as_vector=F, hessian=F, iter=1, verbose=F)
op2.all.pred_delta <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                               intercept=0, beta=op2.all$par$beta_delta),
                             as_vector=F, hessian=F, iter=1, verbose=F)
delta.all.pred.weighted <- stan_matrix_to_datatable(exp(op.all.pred_delta$par$weighted), op.all.pred_delta$par$x)
delta.all.pred <- data.table(pos=op.pred_delta$par$x, delta1=exp(op.all.pred_delta$par$log_mean), delta2=exp(op2.all.pred_delta$par$log_mean))
###plots: op1 vs op2
#nu
a=cbind(biases,data.table(op1=exp(op.all$par$log_nu), 
                          op2=exp(op2.all$par$log_nu), 
                          op.init=exp(op.b1$par$log_nu),
                          delta=op.all$par$log_delta))
ggplot(a[pos>=3166716&pos<=3191637])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op.all$par$eDE+delta)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op.all$par$eDE-delta)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op.all$par$eRJ)),colour="red")+
  geom_point(aes(pos, op1),colour="blue")+
  geom_point(aes(pos, op2),colour="green")+
  geom_line(data=nu.all.pred[pos>=3166716&pos<=3191637], aes(pos, nu1), colour="blue")+
  geom_line(data=nu.all.pred[pos>=3166716&pos<=3191637], aes(pos, nu2), colour="green")+
  geom_line(data=nu.pred[pos>=3166716&pos<=3191637], aes(pos, nu1), colour="yellow")

ggplot(a[pos>=3266716&pos<=3291637])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op.all$par$eDE+delta)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op.all$par$eDE-delta)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op.all$par$eRJ)),colour="red")+
  geom_point(aes(pos, op1),colour="blue")+
  geom_point(aes(pos, op2),colour="green")+
  geom_point(aes(pos, op.init),colour="yellow")+
  geom_line(data=nu.all.pred[pos>=3266716&pos<=3291637], aes(pos, nu1), colour="blue")+
  geom_line(data=nu.all.pred[pos>=3266716&pos<=3291637], aes(pos, nu2), colour="green")+
  geom_line(data=nu.pred[pos>=3266716&pos<=3291637], aes(pos, nu1), colour="yellow")

#delta
a=cbind(biases,data.table(op1=exp(op.all$par$log_delta),
                          op2=exp(op2.all$par$log_delta),
                          log_nu=op.all$par$log_nu))
ggplot(a[pos>=3166716&pos<=3191637])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/(exp(log_nu+op.all$par$eDE))), colour="orange")+
  geom_point(aes(pos, dangling.R/(exp(log_nu+op.all$par$eDE))),  shape=0, colour="pink")+
  geom_point(aes(pos, rejoined/(exp(log_nu+op.all$par$eRJ))), shape=0, colour="red")+
  geom_point(aes(pos, op1),colour="blue")+
  geom_point(aes(pos, op2),colour="green")+
  geom_line(data=delta.all.pred[pos>=3166716&pos<=3191637], aes(pos, delta1), colour="blue")+
  geom_line(data=delta.all.pred[pos>=3166716&pos<=3191637], aes(pos, delta2), colour="green")+
  geom_line(data=delta.pred[pos>=3166716&pos<=3191637], aes(pos, delta1), colour="yellow")
#decay
a=cbind(counts,data.table(decay1=exp(op.all$par$log_decay),
                          decay2=exp(op2.all$par$log_decay),
                          decay.init=exp(op.b3$par$log_decay),
                          cclose1=exp(op.all$par$log_mean_cclose),
                          cclose2=exp(op2.all$par$log_mean_cclose),
                          cfar1=exp(op.all$par$log_mean_cfar),
                          cfar2=exp(op2.all$par$log_mean_cfar),
                          cup1=exp(op.all$par$log_mean_cup),
                          cup2=exp(op2.all$par$log_mean_cup),
                          cdown1=exp(op.all$par$log_mean_cdown),
                          cdown2=exp(op2.all$par$log_mean_cdown)))
ggplot(a[pos2-pos1>1e4][sample(.N,min(.N,10000))])+scale_y_log10()+scale_x_log10()+
  geom_point(aes(pos2-pos1, contact.close*decay1/cclose1), colour="blue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.far*decay1/cclose1), colour="green", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.up*decay1/cclose1), colour="darkblue", alpha=0.01)+
  geom_point(aes(pos2-pos1, contact.down*decay1/cclose1), colour="darkgreen", alpha=0.01)+
  geom_line(aes(pos2-pos1, decay1) ,colour="pink")+
  geom_line(aes(pos2-pos1, decay2) ,colour="orange")+
  geom_line(aes(pos2-pos1, decay.init) ,colour="yellow")







#### subsampling timings
biases=fread("data/caulo_3000000-4000000_biases.dat")
setkey(biases,id)
counts=fread("data/caulo_3000000-4000000_counts.dat")
both=fread("data/caulo_3000000-4000000_both.dat")

system.time(op.full <- optimize_all_noinit(smfit, biases, counts, Krow=1000, Kdiag=10,
                                      lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
system.time(op.100000 <- optimize_all_noinit(smfit, biases, counts[sample(.N,100000)], Krow=1000, Kdiag=10,
                                      lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
system.time(op.60000 <- optimize_all_noinit(smfit, biases, counts[sample(.N,60000)], Krow=1000, Kdiag=10,
                                             lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
system.time(op.30000 <- optimize_all_noinit(smfit, biases, counts[sample(.N,30000)], Krow=1000, Kdiag=10,
                                            lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
system.time(op.10000 <- optimize_all_noinit(smfit, biases, counts[sample(.N,10000)], Krow=1000, Kdiag=10,
                                            lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
system.time(op.5000 <- optimize_all_noinit(smfit, biases, counts[sample(.N,5000)], Krow=1000, Kdiag=10,
                                            lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))
rel_var = function(a,b) {
  return(mean((a-b)/a))
}
max_var = function(a,b) {
  return(max(abs((a-b)/a)))
}
times=data.table(name=c("full","100k", "60k", "30k", "10k", "5k"), nsteps=c(1602, 1779, 1424, 1043, 832, 598),
           time=c(681, 608, 299, 110, 32, 11), npoints=c(counts[,.N], 100000, 60000, 30000, 10000, 5000))
ggplot(times)+geom_point(aes(npoints,time))+scale_x_log10()+scale_y_log10()
lm(data=times, log(time)~log(npoints))

exps=list(n100000=op.100000,n60000=op.60000,n30000=op.30000,n10000=op.10000,n5000=op.5000)
names(op.full$par)
sapply(exps, function(x){rel_var(op.full$par$eC, x$par$eC)})
sapply(exps, function(x){rel_var(op.full$par$eRJ, x$par$eRJ)})
sapply(exps, function(x){mean(op.full$par$log_nu-x$par$log_nu, na.rm=T)})
sapply(exps, function(x){max(abs(op.full$par$log_nu-x$par$log_nu), na.rm=T)})
sapply(exps, function(x){sum(abs(op.full$par$log_nu-x$par$log_nu)>0.01, na.rm=T)/length(op.full$par$log_nu)})\
#output:
#n100000    n60000    n30000    n10000     n5000 
#0.4723810 0.7390476 0.8552381 0.9238095 0.9504762 
sapply(exps, function(x){sum(abs(op.full$par$log_nu-x$par$log_nu)>0.1, na.rm=T)/length(op.full$par$log_nu)})
#output:
#n100000     n60000     n30000     n10000      n5000 
#0.00000000 0.02095238 0.09523810 0.31047619 0.49714286 

a=data.table(sapply(exps, function(x){abs(op.full$par$log_nu-x$par$log_nu)}))
ggplot(melt(a)[value<0.5])+geom_histogram(aes(x=value, y=..density..),binwidth=0.01)+facet_grid(variable ~ .)+
  ggtitle("error on nu by downsampling 125k counts")
ggsave(filename = "error_nu.png", width=5, height=4)
#
a=data.table(sapply(exps, function(x){abs(op.full$par$log_delta-x$par$log_delta)}))
ggplot(melt(a)[value<0.5])+geom_histogram(aes(x=value, y=..density..),binwidth=0.01)+facet_grid(variable ~ .)+
  ggtitle("error on delta by downsampling 125k counts")
ggsave(filename = "error_delta.png", width=5, height=4)
