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
ggplot(matdiff[step>3])+geom_raster(aes(bin1,bin2,fill=-diff))+geom_raster(aes(bin2,bin1,fill=-(diff)))+
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

perf.c = csnorm:::gfl_perf_iteration(ctsg, dispersion, diag.rm, nbins, trails, lambda2, ctsg.ref,
                            alpha=5, inflate=2, tol.value=1e-3, nperf=100, maxsteps=1e5, perf.c)


value.C = csnorm:::gfl_perf_iteration(cts[name!=ref], csg@par$alpha, diag.rm, csg@par$nbins, trails, 5, ref=cts[name==ref],
                                        alpha=5, inflate=2, tol.value=1e-6, nperf=100, maxsteps=100000, state=init.state)
ggplot(as.data.table(value.C$mat)[,.(bin1,bin2,delta=value.C$delta,phi.ref=value.C$phi.ref)])+geom_raster(aes(bin1,bin2,fill=phi.ref))+scale_fill_gradient2()





matg.diff=as.data.table(csnorm:::cts_to_diff_mat(ctsg,ctsg.ref,nbins,dispersion,rep(0,nbins*(nbins+1)/2),rep(0,nbins*(nbins+1)/2),diag.rm))
ggplot(matg.diff)+geom_raster(aes(bin1,bin2,fill=phihat.ref))+geom_raster(aes(bin2,bin1,fill=deltahat))+scale_fill_gradient2()
ggplot(matg.diff)+geom_raster(aes(bin1,bin2,fill=pmin(5,phihat.ref)))+geom_raster(aes(bin2,bin1,fill=pmax(-5,pmin(5,deltahat))))+scale_fill_gradient2()
ggplot(matg.diff)+geom_raster(aes(bin1,bin2,fill=pmax(-5,pmin(5,phihat.ref*weight))))+geom_raster(aes(bin2,bin1,fill=pmax(-5,pmin(5,deltahat*weight))))+scale_fill_gradient2()
matg=as.data.table(csnorm:::cts_to_signal_mat(ctsg,nbins,dispersion,rep(0,nbins*(nbins+1)/2),diag.rm))
ggplot(matg)+geom_raster(aes(bin1,bin2,fill=phihat))+scale_fill_gradient2()
matg.ref=as.data.table(csnorm:::cts_to_signal_mat(ctsg.ref,nbins,dispersion,rep(0,nbins*(nbins+1)/2),diag.rm))
ggplot(matg.ref)+geom_raster(aes(bin1,bin2,fill=phihat))+scale_fill_gradient2()
ggplot()+geom_raster(data=matg,aes(bin1,bin2,fill=pmin(phihat,5)))+geom_raster(data=matg.ref,aes(bin2,bin1,fill=pmin(phihat,5)))+scale_fill_gradient2()
ggplot()+geom_raster(data=matg,aes(bin1,bin2,fill=weight))+geom_raster(data=matg.ref,aes(bin2,bin1,fill=weight))+scale_fill_gradient2()


mat=merge(matg,matg.ref,by=c("bin1","bin2","diag.idx"),suffixes=c("",".ref"))
mat[,c("deltahat","deltahat.var","ncounts.ref","weight.ref"):=list(phihat-phihat.ref,phihat.var+phihat.var.ref,NULL,NULL)]
setcolorder(mat,names(matg.diff))
mat[,c("ncounts","weight"):=NULL]
matg.diff[,c("ncounts","weight"):=NULL]
setkeyv(mat,key(matg.diff))
all.equal(mat,matg.diff)
merge(mat,matg.diff,by=c("bin1","bin2","diag.idx"),suffixes=c(".man",".diff"))[,all.equal(phihat.var.ref.diff,phihat.var.ref.man)]

matg.diff2=as.data.table(csnorm:::cts_to_diff_mat(ctsg.ref,ctsg,nbins,dispersion,rep(0,nbins*(nbins+1)/2),rep(0,nbins*(nbins+1)/2),diag.rm))
setnames(matg.diff2,c("phihat","phihat.var","phihat.ref","phihat.var.ref"),c("phihat.ref","phihat.var.ref","phihat","phihat.var"))
matg.diff2[,deltahat:=-deltahat]
setcolorder(matg.diff2,names(matg.diff))
all.equal(matg.diff,matg.diff2)


matg = csnorm:::gfl_get_matrix(ctsg, nbins, dispersion, diag.rm, trails, 0, lambda2, 0,
                               tol.value=tol.val, state=state, ctsg.ref)
matg2 = csnorm:::gfl_get_matrix(ctsg.ref, nbins, dispersion, diag.rm, trails, 0, lambda2, 0,
                               tol.value=tol.val, state=state, ctsg.ref=ctsg)[,.(bin1,bin2,phihat.ref=)]
ggplot(matg)+geom_raster(aes(bin1,bin2,fill=phihat.ref))+geom_raster(aes(bin2,bin1,fill=valuehat))+scale_fill_gradient2()


#norm
csnorm:::csnorm_fused_lasso(csig, positive=T, fixed=F, constrained=constrained, simplified=T, verbose=verbose)
-> lambda1_eCprime_simplified
#signal
csnorm:::csnorm_fused_lasso(csig, positive=T, fixed=T, constrained=T, simplified=T, verbose=verbose)
-> lambda1_only
#diff
csnorm:::csnorm_fused_lasso(csig, positive=F, fixed=T, constrained=T, simplified=F, verbose=verbose,
                            ctsg.ref=csig@cts.ref)
-> lambda1_only

#negative case
load("tmp_negative_csig.RData")
lambda2 #1.399259
a=csnorm:::gfl_BIC(csig, lambda2)
a$eCprime #0.979505
matg=csnorm:::gfl_get_matrix(csig, a$lambda1, a$lambda2, a$eCprime)
print(ggplot(matg)+geom_raster(aes(bin1,bin2,fill=value))+scale_fill_gradient2())



load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_base10k_bpk30_dfuse20_cpp.RData")
cts = csnorm:::csnorm_gauss_signal_muhat_mean(cs, cs@zeros, cs@settings$sbins)
groupnames=cts[,unique(name)]
g=groupnames[1]
csig=new("CSbsig", mat=cs@par$signal[name==g], trails=cs@settings$trails, cts=cts[name==g],
         settings=list(diag.rm=diag.rm, nbins=nbins, dispersion=cs@par$alpha, tol.val=cs@settings$tol.leg,
                       inflate=2, nperf=500, opt.every=10, maxsteps=100000))
csig@state = csnorm:::gfl_compute_initial_state(csig, diff=F, init.alpha=5)


lambda2=0.5133119
a=csnorm:::gfl_BIC(csig,lambda2=lambda2, constrained=T, positive=T, fixed=F)
a=csnorm:::gfl_get_matrix(csig,0,lambda2,0)
mat=as.data.table(as.data.table(a))[,.(bin1,bin2,phihat=valuehat,ncounts,weight,beta=value,diag.idx)]
mat[,summary(beta)]
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=-beta))+scale_fill_gradient2()
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=weight))

minval=mat[,min(beta)]
test = foreach (eCprime=seq(minval,0,length.out=100), .combine=rbind) %do% {
  lambda1=-eCprime
  mat[,soft:=sign(beta-eCprime)*pmax(0,abs(beta-eCprime)-lambda1)+eCprime]
  mat[,.(eCprime=eCprime, lambda1=lambda1, BIC=sum(weight*(phihat-soft)^2), is.positive=all(soft>=eCprime))]
}
ggplot(test)+geom_line(aes(eCprime,BIC,colour=is.positive))

lvals = foreach(UB=mat[,unique(beta)], .combine=rbind) %dopar% {
  #rm(beta)
  lv=mat[,.(UB=UB,lbound=(UB-min(beta))/2,lmin=weighted.mean(pmax(beta,UB)-phihat, weight),
         lsnc=log(sum(ncounts)), nw=sum(weight))]
  #b=csnorm:::gfl_BIC_fixed(csig, lambda1=lv[,pmax(lmin,lbound)], lambda2=lambda2, eCprime=UB-lv[,pmax(lmin,lbound)])
  b=csnorm:::gfl_BIC_fixed(csig, lambda1=lv[,lbound], lambda2=lambda2, eCprime=UB-lv[,lbound])
  lv[,c("reBIC","relsnc","reR","redof","nub"):=list(b$BIC,log(sum(b$mat$ncounts)),lsnc*b$dof,b$dof,b$dof)]
}
lvals[,c("L","R","Lmod"):=list(nw*(pmax(lmin,lbound)-lmin)^2,lsnc*nub,nw*(lbound-lmin)^2)]
lvals[,c("BIC","BICmod"):=list(L+R,Lmod+R)]
setkey(lvals,UB)
ggplot(melt(lvals[,.(UB,lmin,lbound)],id.vars = "UB"))+geom_line(aes(UB,value,colour=variable))
ggplot(melt(lvals[,.(UB,L,reL=reBIC-reR)],id.vars = "UB"))+geom_line(aes(UB,value,colour=variable))
ggplot(melt(lvals[,.(UB,R,reR)],id.vars = "UB"))+geom_line(aes(UB,value,colour=variable))
ggplot(melt(lvals[,.(UB,nub,redof)],id.vars = "UB"))+geom_line(aes(UB,value,colour=variable))
ggplot(melt(lvals[,.(UB,BIC,reBIC)],id.vars = "UB"))+geom_line(aes(UB,value,colour=variable))#+ylim(0,3000)+geom_line(aes(V5+V7,V9,colour=V1),data=a[V2=="grid"])
ggplot(lvals[,.(UB,BIC,reBIC)])+geom_line(aes(UB,BIC-reBIC))

log(sum(b$mat$ncounts))

lv=mat[,.(UB=UB,lbound=(UB-min(beta))/2,lmin=weighted.mean(pmax(beta,UB)-phihat, weight),
          lsnc=log(sum(ncounts)), nw=sum(weight), nub=uniqueN(c(UB,pmax(beta,UB)))-1)]
lv[,reBIC:=csnorm:::gfl_BIC_fixed(csig, lambda1=lbound, lambda2=lambda2, eCprime=UB-lbound)$BIC]
lv[,c("L","R","Lmod"):=list(nw*(pmax(lmin,lbound)-lmin)^2,lsnc*nub,nw*(lbound-lmin)^2)]
lv[,c("BIC","BICmod"):=list(L+R,Lmod+R)]
lv
a[V2=="grid",.(UB=V5+V7,BIC=V9,ori=V1,lambda1=V5,eCprime=V7)][UB>0.501&UB<0.503]
-> BIC or BICmod are wrong, recalculation provides correct values



a1 = foreach(i=1:a[V2=="grid",.N],.combine=rbind) %dopar% {
  a[V2=="grid"][i,.(ori=V1,status=V3,lambda1=V5,eCprime=V7,BIC=V9,dof=V11,
                    reBIC=csnorm:::gfl_BIC_fixed(csig, lambda1=V5, lambda2=lambda2, eCprime=V7)$BIC)]
}
ggplot(melt(a1[,.(ori,lambda1,eCprime,BIC,reBIC)],id.vars = c("lambda1","eCprime","ori")))+geom_line(aes(lambda1+eCprime,value,colour=variable))

ggplot()+geom_line(data=a[V2=="grid"],aes(V5,V7,colour=V1))


UB=0.398938
#LB=-0.15
a = foreach (LB=seq(-0.5,0,length.out=100),.combine=rbind) %do% {
  as.data.table(csnorm:::gfl_BIC_fixed(csig,lambda1=(UB-LB)/2,lambda2=2,eCprime=(UB+LB)/2)[c("lambda1","lambda2","eCprime","BIC","dof")])
}
ggplot(a)+geom_line(aes(UB-2*lambda1,BIC))


for (resolution in c(5000,10000,20000)) {
  cat("resolution ", resolution, "\n")
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  csi = csnorm:::prepare_signal_estimation(cs, csg, resolution, 1e-4)
  groupnames=csi@cts[,unique(name)][1:2]
  registerDoParallel(cores=4)
  for(g in groupnames) {
    cat("dataset ",as.character(g),"\n")
    csig = csi
    csig@cts = csi@cts[name==g]
    csig@mat = csi@mat[name==g]
    csig@state = csnorm:::gfl_compute_initial_state(csig, diff=F, init.alpha=5)
    obj = function(x) {
        csig@state <<- csnorm:::gfl_BIC(csig, lambda2=10^(x), constrained=constrained, positive=positive, fixed=fixed)
        cat(as.character(g)," grid ", resolution, " lambda2= ",csig@state$lambda2, " lambda1= ",csig@state$lambda1,
          " eCprime= ",csig@state$eCprime," BIC= ",csig@state$BIC, " dof= ",csig@state$dof,"\n")
      return(csig@state$BIC)
    }
    dt.fix = foreach (lam=10^(seq(0,1,length.out=100)),.combine=rbind) %dopar% {obj(log10(lam))}
  }
}

for (resolution in c(5000,10000,20000)) {
  idx1=get_cs_group_idx(cs, resolution, group, raise=T)
  csg=cs@groups[[idx1]]
  csi = csnorm:::prepare_signal_estimation(cs, csg, resolution, 1e-4)
  groupnames=csi@cts[,unique(name)][1:2]
  registerDoParallel(cores=4)
  for(g in groupnames) {
    csig = csi
    csig@cts = csi@cts[name==g]
    csig@mat = csi@mat[name==g]
    csig@state = csnorm:::gfl_compute_initial_state(csig, diff=F, init.alpha=5)
    b = csnorm:::gfl_BIC(csig, lambda2=1, constrained=T, positive=T, fixed=T)
    cat("resolution",resolution,"dataset",as.character(g),"ncounts",sum(b$mat$ncounts),"\n")
  }
}


a=fread("test.dat")
a[,prior:=V13*14.68]
a[,lik:=V11-prior]
ggplot(melt(a[,.(dset=V1,resolution=V3,lambda2=V5,lambda1=V7,prior,lik,BIC=V11)],id.vars=c("dset","resolution","lambda2")))+
  geom_line(aes(lambda2,value,colour=variable))+facet_grid(dset~resolution)

