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








biases=fread("data/rao_HICall_chr19_35000000-36000000_biases.dat")
setkey(biases,id)
counts=fread("data/rao_HICall_chr19_35000000-36000000_counts.dat")

biases=fread("data/caulo_3000000-4000000_biases.dat")
setkey(biases,id)
counts=fread("data/caulo_3000000-4000000_counts.dat")
both=fread("data/caulo_3000000-4000000_both.dat")

biases=fread("data/caulo_all_biases.dat")
setkey(biases,id)
counts=fread("data/caulo_all_counts.dat")
both=fread("data/caulo_all_both.dat")







#make sure the 6-cutter model fits it nicely
sub=both[sample(.N,min(100000,.N))]
fit=gam(N ~ s(log(distance)) + category:(log(dangling.L.1+1) + log(dangling.R.1+1) + log(rejoined.1+1) +
                                           log(dangling.L.2+1) + log(dangling.R.2+1) + log(rejoined.2+1)),
        data=sub, family=nb())
summary(fit)
sub[,fij:=exp(predict.gam(fit, sub, type="terms", terms="s(log(distance))"))]
sub[,mean:=fit$fitted.values]
ggplot(sub[distance>1e4])+geom_point(aes(distance,N*fij/mean),alpha=0.01)+geom_line(aes(distance,fij))+scale_x_log10()+scale_y_log10()



#subsample counts
counts[,prob:=-log(pos2-pos1)]
counts[,prob:=prob/sum(prob)]
counts=counts[sample(.N, min(.N,50000), prob=prob)]



#### simple optimization
system.time(op <- optimize_all_noinit(smfit, biases, counts, Krow=1000, Kdiag=10,
                                      lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))




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
op.b1 <- optimize_nu(smb1, biases, op.b0, Krow=1000, lambda_nu=.1)
op2.b1 <- optimize_nu(smb1, biases, op.b0, Krow=1000, lambda_nu=.1)
#predict nu everywhere
pred = stan_model(file = "predict_spline_sparse.stan")
op.pred_nu <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                        intercept=0, beta=op.b1$par$beta_nu),
                      as_vector=F, hessian=F, iter=1, verbose=F)
op2.pred_nu <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
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
ggplot(a[pos>=3166716&pos<=3191637])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op.b0$par$eDE)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op.b0$par$eDE)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op.b0$par$eRJ)),colour="red")+
  geom_point(aes(pos, op1),colour="blue")+
  geom_point(aes(pos, op2),colour="green")+
  geom_line(data=nu.pred[pos>=3166716&pos<=3191637], aes(pos, nu1), colour="blue")+
  geom_line(data=nu.pred[pos>=3166716&pos<=3191637], aes(pos, nu2), colour="green")+
  geom_line(data=nu.pred.weighted[x>=3166716&x<=3191637], aes(x,value,colour=variable))
message("MAD for nu: ", mad(a$op1-a$op2))
#




#initial guesses for delta, need to set lambda appropriately
smb2 = stan_model(file = "sparse_cs_norm_init_2_delta.stan")
op.b2 <- optimize_delta(smb2, biases, op.b0, op.b1, Krow=1000, lambda_delta=.1)
op2.b2 <- optimize_delta(smb2, biases, op.b0, op.b1, Krow=1000, lambda_delta=.1)
#predict delta everywhere
pred = stan_model(file = "predict_spline_sparse.stan")
op.pred_delta <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                        intercept=0, beta=op.b2$par$beta_delta),
                      as_vector=F, hessian=F, iter=1, verbose=F)
op2.pred_delta <- optimizing(pred, data = list(K=1000, N=10000, xrange=biases[,c(min(pos),max(pos))],
                                         intercept=0, beta=op2.b2$par$beta_delta),
                       as_vector=F, hessian=F, iter=1, verbose=F)
delta.pred.weighted <- stan_matrix_to_datatable(exp(op.pred_delta$par$weighted), op.pred_delta$par$x)
delta.pred <- data.table(pos=op.pred_delta$par$x, delta1=exp(op.pred_delta$par$log_mean), delta2=exp(op2.pred_delta$par$log_mean))
#
a=data.table(name=c("alpha","deviance_proportion_explained"),
             op1=c(op.b2$par$alpha, op.b2$par$deviance_proportion_explained), 
             op2=c(op2.b2$par$alpha, op2.b2$par$deviance_proportion_explained))
#DE op1 vs op2, centered on RJ
a=cbind(biases,data.table(op1=exp(op.b2$par$log_delta),
                          op2=exp(op2.b2$par$log_delta),
                          log_nu=op.b1$par$log_nu))
ggplot(a[pos>=3166716&pos<=3191637])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/(exp(log_nu+op.b0$par$eDE))),colour="orange")+
  geom_point(aes(pos, dangling.R/(exp(log_nu+op.b0$par$eDE))),colour="pink")+
  geom_point(aes(pos, rejoined/(exp(log_nu+op.b0$par$eRJ))),colour="red")+
  geom_point(aes(pos, op1),colour="blue", shape=0)+
  geom_point(aes(pos, op2),colour="green", shape=0)+
  geom_point(aes(pos, 1/op1),colour="blue", shape=0)+
  geom_point(aes(pos, 1/op2),colour="green", shape=0)+
  geom_line(data=delta.pred[pos>=3166716&pos<=3191637], aes(pos, delta1), colour="blue")+
  geom_line(data=delta.pred[pos>=3166716&pos<=3191637], aes(pos, delta2), colour="green")+
  geom_line(data=delta.pred[pos>=3166716&pos<=3191637], aes(pos, 1/delta1), colour="blue", linetype=2)+
  geom_line(data=delta.pred[pos>=3166716&pos<=3191637], aes(pos, 1/delta2), colour="green", linetype=2)+
  geom_line(data=delta.pred.weighted[x>=3166716&x<=3191637], aes(x,value,colour=variable))
#all 3 types
a=cbind(biases,data.table(RJ=exp(op.b0$par$eRJ+op.b1$par$log_nu),
                          DL=exp(op.b0$par$eDE+op.b1$par$log_nu+op.b2$par$log_delta),
                          DR=exp(op.b0$par$eDE+op.b1$par$log_nu-op.b2$par$log_delta)))
a.pred=data.table(pos=nu.pred$pos, RJ=exp(op.b0$par$eRJ)*nu.pred$nu1,
                  DL=exp(op.b0$par$eDE)*nu.pred$nu1*delta.pred$delta1,
                  DR=exp(op.b0$par$eDE)*nu.pred$nu1/delta.pred$delta1)
ggplot(a[pos>=3166716&pos<=3191637])+scale_y_log10()+
  geom_point(aes(pos, rejoined), colour="red")+ geom_point(aes(pos, RJ), shape=0, colour="red")+ geom_line(data=a.pred[pos>=3166716&pos<=3191637], aes(pos, RJ), colour="red")+
  geom_point(aes(pos, dangling.L), colour="orange")+ geom_point(aes(pos, DL), shape=0, colour="orange")+ geom_line(data=a.pred[pos>=3166716&pos<=3191637], aes(pos, DL), colour="orange")+
  geom_point(aes(pos, dangling.R), colour="pink")+ geom_point(aes(pos, DR), shape=0, colour="pink")+ geom_line(data=a.pred[pos>=3166716&pos<=3191637], aes(pos, DR), colour="pink")


#initial guesses for diagonal decay, need to set lambda appropriately
smb3 = stan_model(file = "sparse_cs_norm_init_3_decay.stan")
op.b3 <- optimize_decay(smb3, biases, counts, op.b0, op.b1, op.b2, Kdiag=10, lambda_diag=.1)
op2.b3 <- optimize_decay(smb3, biases, counts, op.b0, op.b1, op.b2, Kdiag=10, lambda_diag=.1)
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
#see how counts correlate with delta
a=melt(cbind(counts[,.(contact.close,contact.far,contact.down,contact.up)],
             data.table(base=op.b3$par$base_count, log_deltai=op.b3$par$log_deltai, log_deltaj=op.b3$par$log_deltaj)),
       measure.vars=c("contact.close","contact.far","contact.down","contact.up"))
fit.ctl = glm(data=a, formula=value ~ 0+variable:(log_deltai+log_deltaj), family=negbin(1))
summary(fit.ctl)
ggplot(a)+geom_point(aes(log_deltai, log(contact.close)))





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











#compare with previous 6-cutter model
fit=gam(N ~ s(log(distance)) + category:(log(dangling.L.1+1) + log(dangling.R.1+1) + log(rejoined.1+1) +
                                         log(dangling.L.2+1) + log(dangling.R.2+1) + log(rejoined.2+1)),
        data=sub, family=nb())
