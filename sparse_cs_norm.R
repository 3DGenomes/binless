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

convert_to_simple = function(data) {
  #convert full cs norm to simple cs norm input
  newdata=data[,.(begin1,begin2,N,rejoined.1,rejoined.2,dangling.L.1,dangling.L.2,dangling.R.1,dangling.R.2,distance)]
  return(newdata[,.(N=sum(N)),by=c("begin1","begin2","rejoined.1","rejoined.2","dangling.L.1","dangling.L.2",
                                   "dangling.R.1","dangling.R.2","distance")])
}

build_stan_list = function(counts, biases, Krow=100, Kdiag=10) {
  return(list(Krow=Krow, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
              danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
              Kdiag=Kdiag, N=counts[,.N],
              counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
              cidx=t(data.matrix(counts[,.(id1,id2)]))))
}

subsample_counts_logspaced = function(counts, sub=1000) {
  prob=counts[,-log(pos2-pos1)]
  prob=prob/sum(prob)
  return(counts[sample(.N, min(.N,sub), prob=prob)])
}

double_genomic_params = function(params, mult=2) {
  params$beta_nu=c(params$beta_nu[1], rep(params$beta_nu, each=mult))
  params$beta_nu_diff=c(params$beta_nu_diff[1:2], rep(params$beta_nu_diff, each=mult))
  params$beta_delta=c(params$beta_delta[1], rep(params$beta_delta, each=mult))
  params$beta_delta_diff=c(params$beta_delta_diff[1:2], rep(params$beta_delta_diff, each=mult))
  return(params)
}

optimized_fit = function(sm, counts, biases, iter=2000, verbose=T, init.counts=10000, mult=2) {
  #fit with a small number of counts first, ~100k resolution
  message("*** Initial fits with increasing basis size and ",init.counts," counts")
  datarange = biases[,max(pos)-min(pos)]
  cts = subsample_counts_logspaced(counts,init.counts)
  data = build_stan_list(cts, biases, Krow=as.integer(datarange/1e5))
  message("*** optimization with Krow=",data$Krow," and N=",data$N)
  op <- optimizing(sm, data = data, as_vector=F, hessian=F, iter=iter, verbose=verbose)
  #increase the number of basis functions by doubling until we reach 1k resolution
  while (datarange/data$Krow>1000) {
    data = build_stan_list(cts, biases, Krow=data$Krow*mult)
    init = double_genomic_params(op$par, mult=mult)
    message("*** optimization with Krow=",data$Krow," and N=",data$N)
    op <- optimizing(sm, data = data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init)
  }
  #do 10, 25 and 100% of counts at 1k res
  message("*** Now increasing the amount of fitted counts")
  for (pc in c(.10,.25,1)) {
    ncounts = as.integer(pc * counts[,.N])
    if (ncounts < cts[,.N]) next
    cts = subsample_counts_logspaced(counts, ncounts)
    data = build_stan_list(cts, biases, Krow=data$Krow)
    message("*** optimization with Krow=",data$Krow," and N=",data$N)
    op <- optimizing(sm, data = data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init)
  }
  return(op)
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

counts[,prob:=-log(pos2-pos1)]
counts[,prob:=prob/sum(prob)]
counts=counts[sample(.N, min(.N,50000), prob=prob)]

#fit it with stan and gam
sm = stan_model(file = "sparse_cs_norm_fit.stan")

#iterative optimization
counts[,dbin:=cut(log(pos2-pos1),10, ordered_result = T, right=F)]
ggplot(subsample_counts_logspaced(counts,50000)[,.N,by=dbin])+geom_bar(aes(dbin,log(N)), stat = "identity")
system.time(op <- optimized_fit(sm, counts, biases))


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
a=cbind(biases,data.table(op1=exp(op.all$par$log_nu), op2=exp(op2.all$par$log_nu)))
ggplot(a[pos>=3166716&pos<=3191637])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(op.all$par$eDE)),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(op.all$par$eDE)),colour="pink")+
  geom_point(aes(pos, rejoined/exp(op.all$par$eRJ)),colour="red")+
  geom_point(aes(pos, op1),colour="blue")+
  geom_point(aes(pos, op2),colour="green")+
  geom_line(data=nu.all.pred[pos>=3166716&pos<=3191637], aes(pos, nu1), colour="blue")+
  geom_line(data=nu.all.pred[pos>=3166716&pos<=3191637], aes(pos, nu2), colour="green")+
  geom_line(data=nu.pred[pos>=3166716&pos<=3191637], aes(pos, nu1), colour="red")
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










#simple optimization
system.time(op <- optimize_all_noinit(smfit, biases, counts, Krow=1000, Kdiag=10,
                                      lambda_nu=.2, lambda_delta=.2, lambda_diag=.1))

#predict fits
smp = stan_model(file="sparse_cs_norm_predict.stan")
system.time(pred <- optimizing(smp, data = list(Krow=4000, Kdiag=10, cutsites=biases[,c(min(pos),max(pos))],
                                                N=10000, S=100000, intercept=op$par$intercept,
                                                eRJ=op$par$eRJ, eDE=op$par$eDE, beta_nu=op$par$beta_nu,
                                                beta_delta=op$par$beta_delta, beta_diag=op$par$beta_diag),
                             as_vector=F, hessian=F, iter=1, verbose=T))
biases.pred=data.table(pos=pred$par$genome, nu=exp(pred$par$log_nu), delta=exp(pred$par$log_delta),
                       mean_DL=pred$par$mean_DL, mean_DR=pred$par$mean_DR, mean_RJ=pred$par$mean_RJ)
counts.pred=data.table(distance=pred$par$distnce, decay=exp(pred$par$log_decay))

#compare with previous 6-cutter model
fit=gam(N ~ s(log(distance)) + category:(log(dangling.L.1+1) + log(dangling.R.1+1) + log(rejoined.1+1) +
                                         log(dangling.L.2+1) + log(dangling.R.2+1) + log(rejoined.2+1)),
        data=sub, family=nb())
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



ggplot(counts[pos2-pos1>1e4][sample(.N,100000)])+scale_y_log10() +scale_x_log10() +
  geom_point(aes(pos2-pos1,contact.up/mean_cup1*decay1), alpha=0.01) +
  geom_line(aes(pos2-pos1,decay1), colour="red") + geom_line(aes(pos2-pos1,decay.gam),colour="green")
ggplot(counts)+scale_y_log10() + geom_point(aes(pos2-pos1,contact.up/mean_cup.gam*decay.gam), alpha=0.1) +
  geom_line(aes(pos2-pos1,decay1), colour="red") + geom_line(aes(pos2-pos1,decay.gam),colour="green")

ggplot(biases) + scale_y_log10() +
  geom_point(aes(pos,dangling.L),colour="red") +
  geom_point(aes(pos,dangling.R),colour="green") +
  geom_point(aes(pos,rejoined),colour="blue") +
  geom_line(data=biases.pred[mean_DL<biases[,max(dangling.L)]&mean_DL>.1], aes(pos,mean_DL),colour="red")+
  geom_line(data=biases.pred[mean_DR<biases[,max(dangling.R)]&mean_DR>.1], aes(pos,mean_DR),colour="green")+
  geom_line(data=biases.pred[mean_RJ<biases[,max(rejoined)]&mean_RJ>.1], aes(pos,mean_RJ),colour="blue")+xlim(3500000,3550000)

ggplot(melt(biases[,.(id,pos,GC.R.400,GC.L.400,GC.R.10,GC.L.10)], id.vars = c("id","pos"))) +
  geom_point(aes(pos, value, colour=variable)) +
  geom_line(aes(pos, value, colour=variable))
pairs(~log(1+dangling.L)+log(1+dangling.R)+log(1+rejoined)+GC.R.400+GC.L.400+GC.R.10+GC.L.10, data=biases)               
  
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



