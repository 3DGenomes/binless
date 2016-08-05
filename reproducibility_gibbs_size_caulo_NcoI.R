library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)
library(Hmisc)
library(dplyr)
library(rstan)

setwd("/home/yannick/simulations/cs_norm")

#zoom on a portion of the dataset
load("data/caulo_NcoI_all_csdata_with_data.RData")
cs.ori=cs

registerDoParallel(cores=30)
foreach (i=seq(1,4,by=0.5)) %dopar% {
  data=cs.ori@data[re.closest1<=i*1000000&re.closest2<=i*1000000]
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=4042929)
  csd = new("CSdata", info=cs@info, settings=cs@settings,
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  cs = merge_cs_norm_datasets(list(csd))
  save(cs, file=paste0("data/caulo_NcoI_",i,"M_csnorm.RData"))
}


### reproducibility: initial guess vs full estimation
sizes=seq(1,4,by=0.5)
lambdas=c(0.01,0.05,0.1,0.5,1)

registerDoParallel(cores=30)
foreach (size=sizes) %:% foreach (lambda=lambdas) %dopar% {
  load(paste0("data/caulo_NcoI_",size,"M_csnorm.RData"))
  cs = run_gibbs(cs, design=NULL, bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10, lambda=lambda,
                       ngibbs = 1, iter=100000)
  save(cs, file=paste0("data/caulo_NcoI_",size,"M_lambda",lambda,"_gibbs1_csnorm_optimized.RData"))
}



# generate plots
outputs = foreach (i=sizes, .combine=rbind) %:% foreach (j=lambdas, .combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"M_lambda",j,"_gibbs1_csnorm_optimized.RData"))
  data.table(dset=paste("sizes",i,"lambda",j),size=i, lambda=j, logp=cs@par$value,
             runtime=cs@par$runtime,
             lambda_nu=cs@par$lambda_nu, lambda_delta=cs@par$lambda_delta, out=tail(cs@par$output,1))
}
outputs
summary(lm(data=outputs, log(runtime)~log(size)+log(lambda)))
ggplot(outputs)+geom_point(aes(size,runtime,colour=log(lambda)))+scale_x_log10()+scale_y_log10()+
  ylab("run time (s)")+xlab("dataset size")+ggtitle("t ~ size^-0.36")

#nu and delta: plots and correlation
registerDoParallel(cores=30)
nu = foreach (i=sizes, .combine=rbind) %:% foreach (j=lambdas, .combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"M_lambda",j,"_gibbs1_csnorm_optimized.RData"))
  rbind(csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$init$beta_nu, beta_delta=cs@par$init$beta_delta,
                                   bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                     ,.(pos,log_nu,log_delta,size=i, lambda=j, status="init")],
        csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                         bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                           ,.(pos,log_nu,log_delta,size=i, lambda=j, status=strsplit(tail(cs@par$output,1)," ")[[1]][3])])
}
ggplot(nu)+scale_y_log10(limits=c(0.1,10))+geom_line(aes(pos,exp(log_nu),colour=status))+facet_grid(size~lambda)+ylab("nu")
ggplot(nu)+scale_y_log10(limits=c(0.1,10))+geom_line(aes(pos,exp(log_delta),colour=status))+facet_grid(size~lambda)+ylab("delta")
ggplot(nu[pos>150000&pos<300000])+scale_y_log10()+geom_line(aes(pos,exp(log_nu),colour=status))+facet_grid(size~lambda)+ylab("nu")
ggplot(nu[pos>150000&pos<300000])+scale_y_log10()+geom_line(aes(pos,exp(log_delta),colour=status))+facet_grid(size~lambda)+ylab("delta")


#diagonal decay: simplified and full compared, all lambdas
registerDoParallel(cores=30)
decay=rbind(foreach (i=sizes,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
             size=i,lambda=j,key="dist")
},
foreach (i=sizes,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
  data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
             size=i,lambda=j,key="dist")
})
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_grid(lambda~size, labeller=label_both)+scale_y_log10()+scale_x_log10()
ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_decayall.png", width=10, height=7.5)

#diagonal decay: optimal lambdas
registerDoParallel(cores=30)
dsimplified = foreach (i=sizes,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized_init.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="simplified",
               size=i,lambda=j,key="dist")
  }
dfull = foreach (i=sizes,.combine=rbind) %:%
  foreach (j=lambdas,.combine=function(x,y){if (x$logp[1]<y$logp[1]){return(y)}else{return(x)}}) %dopar% {
    load(paste0("data/caulo_NcoI_150k_sub",i,"k_lambda",j,"_csnorm_optimized.RData"))
    data.table(dist=cs@counts[,distance],decay=exp(cs@par$log_decay),logp=cs@par$value,cat="full",
               size=i,lambda=j,key="dist")
  }
decay=rbind(dsimplified,dfull)
ggplot(decay)+geom_line(aes(dist,decay,colour=cat))+facet_wrap(~size, labeller=label_both)+
  ylab("decay")+scale_y_log10()+scale_x_log10()
ggsave(filename="images/caulo_NcoI_150k_full_vs_gibbs_decay.png", width=10, height=7.5)



#parameters
registerDoParallel(cores=30)
params = foreach (i=sizes,.combine=rbind) %:% foreach (j=lambdas,.combine=rbind) %dopar% {
  load(paste0("data/caulo_NcoI_",i,"M_lambda",j,"_gibbs1_csnorm_optimized.RData"))
  rbind(data.table(logp=NA,cat="init",size=i,lambda=j,status=NA,
             eC=cs@par$init$eC, eRJ=cs@par$init$eRJ, eDE=cs@par$init$eDE, lambda_diag=cs@par$init$lambda_diag,
             lambda_delta=cs@par$init$lambda_delta, lambda_nu=cs@par$init$lambda_nu, alpha=cs@par$init$alpha),
        data.table(logp=cs@par$value,cat="gibbs",size=i,lambda=j,status=strsplit(tail(cs@par$output,1)," ")[[1]][3],
               eC=cs@par$eC, eRJ=cs@par$eRJ, eDE=cs@par$eDE, lambda_diag=cs@par$lambda_diag,
               lambda_delta=cs@par$lambda_delta, lambda_nu=cs@par$lambda_nu, alpha=cs@par$alpha))
}
setkey(params, size,lambda,cat)
ggplot(params,aes(size,eC,colour=cat))+geom_point()+ylim(-5,20)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_jitter()#+scale_x_log10()
ggplot(params)+geom_point(aes(factor(paste(size,lambda)),lambda_nu,colour=cat))+scale_y_log10()
ggplot(params)+geom_point(aes(factor(paste(size,lambda)),lambda_delta,colour=cat))#+ylim(0,25)
ggplot(dcast(params[,.(cat,size,lambda,eC)],size+lambda~cat)[,.(size,lambda,ratio=init/gibbs)][order(ratio)])+
  geom_point(aes(size,ratio,colour=log(lambda)))+ylim(0,2)
dcast(params[,.(cat,size,lambda,log(lambda_nu))],size+lambda~cat)



size=4
load(paste0("data/caulo_NcoI_",size,"M_csnorm.RData"))
cs@settings = c(cs@settings, list(bf_per_kb=0.25, bf_per_decade=5, bins_per_bf=10, groups=10,
                                  lambda=0.1, iter=1000000))
cs@counts = fill_zeros(counts = cs@counts, biases = cs@biases)
cs@counts[,distance:=pmin(abs(pos2-pos1), cs@settings$circularize+1-abs(pos2-pos1))]
dmin=0.99
dmax=cs@settings$circularize/2+0.01
cs@settings$dmin=dmin
cs@settings$dmax=dmax
#init.a=system.time(init.output <- capture.output(init.par <- run_split_parallel_initial_guess(
#  counts=cs@counts, biases=cs@biases, lambda=0.1, verbose=T,
#  bf_per_kb=0.25, dmin=dmin, dmax=dmax, bf_per_decade=5, iter=1000000)))
init.a=system.time(init.output <- capture.output(init.op <- csnorm:::csnorm_simplified_guess(
  biases = cs@biases, counts = cs@counts, dmin = cs@settings$dmin, dmax = cs@settings$dmax, lambda=cs@settings$lambda,
  groups = cs@settings$groups, bf_per_kb = cs@settings$bf_per_kb, bf_per_decade=cs@settings$bf_per_decade,
  iter = cs@settings$iter)))
init.par=init.op$par
a=system.time(output <- capture.output(op.gen <- csnorm:::csnorm_simplified_genomic(
  biases = cs@biases, counts = cs@counts,
  log_decay = init.par$log_decay, log_nu = init.par$log_nu, log_delta = init.par$log_delta,
  groups = cs@settings$groups, bf_per_kb = cs@settings$bf_per_kb, iter = cs@settings$iter, init=init.par,
  verbose=T, sample_file = paste0("tmp/test_",size,"M_lambda",cs@settings$lambda,"_gibbs.dat"))))
cs@par=op.gen$par
nu=rbind(csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=init.par$beta_nu, beta_delta=init.par$beta_delta,
                                       bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                         ,.(pos,log_nu,log_delta, size=size, lambda=cs@settings$lambda, stage="init",
                                            status=strsplit(tail(init.output,1)," ")[[1]][3])],
      csnorm:::generate_genomic_biases(biases=cs@biases, beta_nu=cs@par$beta_nu, beta_delta=cs@par$beta_delta,
                                       bf_per_kb=cs@settings$bf_per_kb, points_per_kb = 10)[
                                         ,.(pos,log_nu,log_delta,size=size, lambda=cs@settings$lambda, stage="gibbs",
                                            status=strsplit(tail(output,1)," ")[[1]][3])])
ggplot(nu)+scale_y_log10(limits=c(0.1,10))+geom_line(aes(pos,exp(log_nu),colour=stage,linetype=status))+facet_grid(size~lambda)+ylab("nu")
ggplot(nu[pos>150000&pos<300000])+scale_y_log10(limits=c(0.1,10))+geom_line(aes(pos,exp(log_nu),colour=stage,linetype=status))+facet_grid(size~lambda)+ylab("nu")

#look at a few parameters
data=fread("tmp/test_4M_lambda0.1_gibbs.dat", skip=16)
nms=c("lambda_nu", "lambda_delta", "alpha", "eC", "eRJ", "eDE",
      "log_nu.1", "log_delta.1", "log_nu.10", "log_delta.10", "log_nu.100", "log_delta.100",
      "beta_nu.1", "beta_delta.1", "beta_nu.10", "beta_delta.10", "beta_nu.100", "beta_delta.100")
mdata=data[,mget(nms)]
mdata[,time:=.I]
mdata=melt(mdata,id.vars="time")
ggplot(mdata)+geom_point(aes(time,value))+geom_line(aes(time,value))+facet_wrap(~variable,scales="free_y")


#manual optimization for genomic
csub=copy(cs@counts)
csub[,decay:=exp(init.par$log_decay)]
bsub=cs@biases[,.(id)]
bsub[,c("nu","delta"):=list(exp(init.par$log_nu),exp(init.par$log_delta))]
csub=merge(bsub[,.(id1=id,nu,delta)],csub,by="id1",all.x=F,all.y=T)
csub=merge(bsub[,.(id2=id,nu,delta)],csub,by="id2",all.x=F,all.y=T, suffixes=c("2","1"))
#collect all counts on left/right side and put into quantile groups
csb=rbind(csub[,.(pos=pos1,ldist=log(distance),R=(contact.close+contact.down),L=(contact.far+contact.up),others=decay*nu2*(delta2+1/delta2))],
         csub[,.(pos=pos2,ldist=log(distance),R=(contact.far+contact.down),L=(contact.close+contact.up),others=decay*nu1*(delta1+1/delta1))])
setkey(csb,pos)
csb[,cbin:=ntile(L+R,10),by=pos]
csl=dcast(csb[,.(pos,cbin,L)], pos~cbin, value.var="L", fun.aggregate=sum)
csl[,pos:=NULL]
csr=dcast(csb[,.(pos,cbin,R)], pos~cbin, value.var="R", fun.aggregate=sum)
csr[,pos:=NULL]
cso=dcast(csb[,.(pos,cbin,others)], pos~cbin, value.var="others", fun.aggregate=sum)
cso[,pos:=NULL]
#run optimization
Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
data=list(Krow=Krow, S=cs@biases[,.N],
          cutsites=cs@biases[,pos], rejoined=cs@biases[,rejoined],
          danglingL=cs@biases[,dangling.L], danglingR=cs@biases[,dangling.R],
          G=10, counts_sum_left=csl, counts_sum_right=csr, log_decay_sum=log(cso))
save(data, init.par, file = paste0("tmp/test_",size,"M_lambda",cs@settings$lambda,"_gibbs_data.dat"))
op=optimizing(csnorm:::stanmodels$simplified_genomic, data=data, as_vector=F, hessian=F,
           iter=100, verbose=T, init=init.par, algorithm="LBFGS", init_alpha=1e-5)



