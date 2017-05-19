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
niter=10
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


value.c = csnorm:::wgfl_perf(cts, cs@par$alpha, niter, trails$nrow, trails$ntrails, trails$trails,
                             trails$breakpoints, lambda2, alpha, inflate, maxsteps, tol.val/2)

all.equal(value.c,value.r)



### plots

plot_diagnostics(cs)$plot
plot_diagnostics(cs)$plot2

signals=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  if ("signal" %in% cs@diagnostics$params[step==i,leg]) {
    sig=copy(cs@diagnostics$params[step==i&leg=="signal",signal][[1]])
    sig[,step:=i]
    sig
  }
}

ggplot(signals[name==name[1]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")
ggplot(signals[step>=step[.N]-1])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_grid(step~ name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")
ggplot(signals[step>=step[.N]])+geom_raster(aes(bin1,bin2,fill=phi))+geom_raster(aes(bin2,bin1,fill=phi))+facet_wrap(~name)+scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")


matdiff=merge(signals[step>step[1],.(name,step,bin1,bin2,phi)],
              signals[step<step[.N],.(name,step=step+1,bin1,bin2,phi)],
              by=c("name","step","bin1","bin2"))[,.(name,step,bin1,bin2,diff=phi.x-phi.y)]
ggplot(matdiff[step>30])+geom_raster(aes(bin1,bin2,fill=-diff))+geom_raster(aes(bin2,bin1,fill=-(diff)))+
  scale_fill_gradient2()+facet_wrap(~step)+coord_fixed()



decays=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  sig=copy(cs@diagnostics$params[step==i&leg=="decay",decay][[1]])
  sig[,step:=i]
  sig
}
ggplot(decays[name==name[1]])+geom_pointrange(aes(distance,kappahat,ymin=kappahat-std,ymax=kappahat+std),alpha=0.1)+
  geom_line(aes(distance,kappa))+facet_wrap(~ step)+scale_x_log10()#+
#geom_line(aes(distance,base_count),colour="red",data=cs@counts[name==name[1]][sample(.N,min(.N,10000))])
ggplot(decays[name==name[1]&step>20])+geom_line(aes(distance,kappa,group=step,colour=factor(step)))+scale_x_log10()

decaydiff=merge(decays[step>step[1],.(name,step,distance,log_decay)],
              decays[step<step[.N],.(name,step=step+1,distance,log_decay)],
              by=c("name","step","distance"))[,.(name,step,distance,diff=log_decay.x-log_decay.y)]
ggplot(decaydiff[step>20])+geom_line(aes(distance,diff,group=step,colour=factor(step)))+scale_x_log10()+facet_wrap(~step)



load("tmp_oscillations.RData")
cs = run_gauss(cs, restart=F, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 10, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-3, tol.leg=1e-3,
               ncounts = 1000000, ncores=ncores, base.res=5000, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)



#works
cs@settings$tol.leg=1e-5
cs@settings$tol.obj=1e-3
cs = run_gauss(cs, restart=T, bf_per_kb=bpk, bf_per_decade=bpd, bins_per_bf=bpb,
               ngibbs = 10, iter=100000, init_alpha=1e-7, init.dispersion = 1, tol.obj=1e-1, tol.leg=1e-3,
               ncounts = 1000000, ncores=ncores, base.res=5000, fit.signal=T, fit.disp=T, fit.decay=T, fit.genomic=T)

#doesnt
cs@settings$tol.decay=1e-3
cs@settings$tol.signal=1e-5
cs@settings$tol.val=1e-1
cs@settings$tol.obj=1e-1


mat=merge(mat.decay,mat.nodecay,suffixes=c(".d",".nd"),by=c("name","bin1","bin2","diag.idx"))
ggplot(mat)+facet_wrap(~name)+geom_raster(aes(bin1,bin2,fill=value.d-value.nd))+geom_raster(aes(bin2,bin1,fill=value.d-value.nd))+
  scale_fill_gradient2()
ggplot(mat)+geom_point(aes(weight.d,weight.nd,colour=weight.d!=weight.nd))+stat_function(fun=identity)

nbins=length(cs@settings$sbins)-1
dispersion=cs@par$alpha
trails=cs@settings$trails
tol.val=cs@settings$tol.signal
csnorm:::optimize_lambda2(cts.nodecay, nbins, dispersion, trails, tol.val = tol.val)



sub="SELP_150k"
bpk=30
dfuse=20
load(paste0("data/rao_HiCall_",sub,"_csnorm_optimized_base5k_bpk",bpk,"_dfuse",dfuse,".RData"))
resolution=5000
group="all"
eCprime=cs@design[,.(eCprime=0),keyby=name]
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores, niter=1)
idx1=get_cs_group_idx(cs, resolution, group, raise=T)
csg=cs@groups[[idx1]]
diag.rm = ceiling(csg@par$dmin/resolution)
cts=csg@cts[,.(name,bin1,bin2,count,lmu.nosig,weight,phi,var)]
ref=cs@design[.N,name]

###signal
stuff = csnorm:::prepare_signal_matrix(cs, csg@names, resolution)
mat=stuff$mat
trails=stuff$trails
#R signal
mat.R = csnorm:::csnorm_compute_raw_signal(cts[name==ref], csg@par$alpha, mat[name==ref], eCprime)
mat.R = mat.R[,.(bin1,bin2,phihat,phihat.var,ncounts,weight,diag.idx)]
mat.R[diag.idx<=diag.rm,c("phihat.var","weight"):=list(Inf,0)]
#C signal
mat.C = as.data.table(csnorm:::cts_to_signal_mat(cts[name==ref], csg@par$nbins, csg@par$alpha,
                                               mat[name==ref,phi], diag.rm))
all.equal(mat.C,mat.R)
merge(mat.C,mat.R,by=c("bin1","bin2"),suffixes=c(".R",".C"))[,summary(phihat.R-phihat.C)]
merge(mat.C,mat.R,by=c("bin1","bin2"),suffixes=c(".R",".C"))[phihat.R-phihat.C==-1]#[,.(unclass(bin1),unclass(bin2))]


###diff
stuff = csnorm:::prepare_difference_matrix(cs, csg@names, resolution, ref)
mat=stuff$mat
trails=stuff$trails
#R diff
mat.R = csnorm:::csnorm_compute_raw_differential(csg@cts, csg@par$alpha, mat, eCprime, ref)
mat.R = mat.R[,.(bin1,bin2,phihat,phihat.var,phihat.ref,phihat.var.ref,deltahat,deltahat.var,ncounts,weight,diag.idx)]
#C diff
mat.C = as.data.table(csnorm:::cts_to_diff_mat(csg@cts[name!=ref], csg@cts[name==ref], csg@par$nbins, csg@par$alpha,
  mat[name!=ref,phi.ref], mat[name!=ref,delta], diag.rm))

all.equal(mat.C,mat.R)
merge(mat.C,mat.R,by=c("bin1","bin2"),suffixes=c(".R",".C"))[,summary(deltahat.R-deltahat.C)]


value.C = csnorm:::gfl_perf_iteration(cts[name!=ref], csg@par$alpha, diag.rm, csg@par$nbins, trails, 5, ref=cts[name==ref],
                                        alpha=5, inflate=2, tol.value=1e-6, nperf=100, maxsteps=100000, state=NULL)
ggplot(as.data.table(value.C$mat)[,.(bin1,bin2,delta=value.C$delta,phi.ref=value.C$phi.ref)])+geom_raster(aes(bin1,bin2,fill=phi.ref))+scale_fill_gradient2()

