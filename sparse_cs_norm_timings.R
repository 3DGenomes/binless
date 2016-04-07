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

counts[,prob:=-log(pos2-pos1)]
counts[,prob:=prob/sum(prob)]
counts=counts[sample(.N, min(.N,50000), prob=prob)]

#fit it with stan and gam
sm = stan_model(file = "sparse_cs_norm_fit.stan")

system.time(op <- optimizing(sm, data = list(Krow=500, S=biases[,.N], cutsites=biases[,pos], rejoined=biases[,rejoined],
                                             danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
                                             Kdiag=10, N=counts[,.N],
                                             counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
                                             cidx=t(data.matrix(counts[,.(id1,id2)]))),
                             as_vector=F, hessian=F, iter=2000, verbose=T, init_alpha=1))

timings.N=data.table(N=c(124000,100000,80000,60000,40000,20000,10000,50000),tnol=c(110,65,65,61,24,15,5,34), t=c(106,103,60,46,35,33,12,43))
summary(lm(data=timings.N, log(tnol)~log(N)))
ggplot(timings.N,aes(N,tnol))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+ggtitle("lambda_bias fixed: t ~ N^1.17")
ggsave(filename = "images/caulo_3000000-4000000_timing_N_lambda_fixed.png", width=10, height=7.5)
summary(lm(data=timings.N, log(t)~log(N)))
ggplot(timings.N,aes(N,t))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+ggtitle("lambda_bias fixed: t ~ N^0.8")
ggsave(filename = "images/caulo_3000000-4000000_timing_N.png", width=10, height=7.5)

#
timings.Krow=data.table(Krow=c(100,50,25,150,200,250,400,500,750,1000), tnol=c(34,33,34,38,46,49,56,51,69,89), t=c(43,38,37,38,38,57,191,252,252,252))
summary(lm(data=timings.Krow, log(t)~log(Krow)))
ggplot(timings.Krow,aes(Krow,t))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+ggtitle("lambda_bias fixed: t ~ Krow^0.25")
ggsave(filename = "images/caulo_3000000-4000000_timing_Krow_lambda_fixed.png", width=10, height=7.5)
summary(lm(data=timings.Krow[1:6], log(t)~log(Krow)))
ggplot(timings.Krow,aes(Krow,t))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+ggtitle("lambda_bias fixed: t ~ Krow^0.25")
ggsave(filename = "images/caulo_3000000-4000000_timing_Krow.png", width=10, height=7.5)


