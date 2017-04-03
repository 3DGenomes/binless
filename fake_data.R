library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)

setwd("/home/yannick/simulations/cs_norm")

csd=generate_fake_dataset(signal=F,eC=-4.4,replicate="1",condition="WT")
save(csd,file="data/fake_replicate1_csdata.RData")
biases.ref=csd@biases
csd=generate_fake_dataset(signal=F,biases.ref=biases.ref,eC=-5,replicate="2",condition="WT")
save(csd,file="data/fake_replicate2_csdata.RData")
csd=generate_fake_dataset(signal=F,eC=-4,biases.ref=biases.ref,replicate="3",condition="WT")
save(csd,file="data/fake_replicate3_csdata.RData")
csd=generate_fake_dataset(signal=F,eC=-3.8,biases.ref=biases.ref,replicate="4",condition="WT")
save(csd,file="data/fake_replicate4_csdata.RData")
csd=generate_fake_dataset(signal=F,eC=-4.1,biases.ref=biases.ref,replicate="5",condition="WT")
save(csd,file="data/fake_replicate5_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.3,replicate="1",condition="KO")
save(csd,file="data/fake_signal_replicate1_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-6.1,replicate="2",condition="KO")
save(csd,file="data/fake_signal_replicate2_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-3.7,replicate="3",condition="KO")
save(csd,file="data/fake_signal_replicate3_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-4.2,replicate="4",condition="KO")
save(csd,file="data/fake_signal_replicate4_csdata.RData")
csd=generate_fake_dataset(signal=T,biases.ref=biases.ref,eC=-5.2,replicate="5",condition="KO")
save(csd,file="data/fake_signal_replicate5_csdata.RData")

csd@counts[,bin1:=round(pos1/10000)]
csd@counts[,bin2:=round(pos2/10000)]
binned=csd@counts[,.(count=sum(contact.far+contact.close+contact.up+contact.down)),by=c("bin1","bin2")]
ggplot(binned)+geom_raster(aes(bin1,bin2,fill=log(count)))

load("data/fake_csnorm_optimized.RData")

load("data/fake_replicate1_csdata.RData")
csd1=csd
load("data/fake_replicate2_csdata.RData")
csd2=csd
load("data/fake_replicate3_csdata.RData")
csd3=csd
load("data/fake_replicate4_csdata.RData")
csd4=csd
load("data/fake_replicate5_csdata.RData")
csd5=csd
cs=merge_cs_norm_datasets(list(csd1,csd2,csd3,csd4,csd5), different.decays="none")
#cs=merge_cs_norm_datasets(list(csd1), different.decays="none")
cs = run_gauss(cs, restart=F, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10, ngibbs = 5, base.res=20000,
               iter=100000, init_alpha=1e-7, ncounts = 100000, type="perf", ncores=30, fit.signal=F)
cs = run_gauss(cs, restart=T, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10, ngibbs = 1, base.res=20000,
               iter=100000, init_alpha=1e-7, ncounts = 100000, type="perf", ncores=30, fit.signal=T)
#cs@par$signal[,phi:=2*phi]
cs = run_gauss(cs, restart=T, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10, ngibbs = 5, base.res=20000,
               iter=100000, init_alpha=1e-7, ncounts = 100000, type="perf", ncores=30, fit.signal=F)
#save(cs,file="data/fake_signal_shrink10pc_new_csnorm_optimized.RData")
cs = run_gauss(cs, restart=T, bf_per_kb=30, bf_per_decade=10, bins_per_bf=10, ngibbs = 3, base.res=20000,
               iter=100000, init_alpha=1e-7, ncounts = 100000, type="perf", ncores=30, fit.signal=T)
#save(cs,file="data/fake_replicate1_signal_shrink10pc_new2_csnorm_optimized.RData")

save(cs,file="data/fake_csnorm_optimized.RData")


plot_diagnostics(cs)$plot
plot_diagnostics(cs)$plot2

ggplot(cs@par$decay)+geom_pointrange(aes(distance,kappahat,ymin=kappahat-std,ymax=kappahat+std),alpha=0.1)+
  geom_line(aes(distance,kappa))+facet_wrap(~ name)+scale_x_log10()+
  geom_line(aes(distance,base_count),colour="red",data=cs@counts[sample(.N,min(.N,100000))])
ggplot(cs@par$decay[distance<1e4])+geom_pointrange(aes(distance,kappahat,ymin=kappahat-std,ymax=kappahat+std),alpha=0.1)+
  geom_line(aes(distance,kappa))+facet_wrap(~ name)+scale_x_log10()+
  geom_line(aes(distance,base_count),colour="red",data=cs@counts[sample(.N,min(.N,100000))][distance<1e4])

ggplot(cs@par$biases[cat=="dangling L"])+geom_point(aes(pos,etahat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_wrap(~ name)+xlim(550000,650000)+
  geom_line(aes(pos,true_log_mean_DL),colour="red",data=cs@biases)

ggplot(cs@par$biases[cat=="dangling R"])+geom_point(aes(pos,etahat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_wrap(~ name)+xlim(550000,650000)+
  geom_line(aes(pos,true_log_mean_DR),colour="red",data=cs@biases)

true_sig=cs@counts[,.(name,pos1,pos2,true_phi)]
true_sig[,bin1:=cut(pos1,cs@settings$sbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
true_sig[,bin2:=cut(pos2,cs@settings$sbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
true_sig=true_sig[,.(phi=mean(true_phi)),keyby=c("name","bin1","bin2")]
ggplot(true_sig[name==name[1]])+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
cs@par$signal=true_sig[cs@par$signal[,.(name,bin1,bin2)]]

#signals
signals=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  if ("signal" %in% cs@diagnostics$params[step==i,leg]) {
    sig=copy(cs@diagnostics$params[step==i&leg=="signal",signal][[1]])
    sig[,step:=i]
    sig
  }
}
ggplot(signals[name==name[1]])+geom_raster(aes(bin1,bin2,fill=phi))+facet_wrap(~ step)+scale_fill_gradient2()

#decays
decays=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  sig=copy(cs@diagnostics$params[step==i&leg=="decay",decay][[1]])
  sig[,step:=i]
  sig
}
ggplot(decays[name==name[1]])+geom_pointrange(aes(distance,kappahat,ymin=kappahat-std,ymax=kappahat+std),alpha=0.01)+
  geom_line(aes(distance,kappa))+facet_wrap(~ step)+scale_x_log10()+
  geom_line(aes(distance,base_count),colour="red",data=cs@counts[name==name[1]][sample(.N,min(.N,10000))])

#biases
biases=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  sig=copy(cs@diagnostics$params[step==i&leg=="bias",biases][[1]])[cat%in%c("contact L", "contact R")]
  delta=c(-4.3,-6.3,-3.7,-3.1,-4.2)-1#cs@diagnostics$params[step==i&leg=="bias",eC[[1]]-eDE[[1]]]
  sig[,step:=i]
  tcl=merge(cs@biases,cbind(cs@design[,.(name)],delta=delta),by="name")[,.(name,id,cat="contact L",true_eta=true_log_mean_DL+delta)]
  tcr=merge(cs@biases,cbind(cs@design[,.(name)],delta=delta),by="name")[,.(name,id,cat="contact R",true_eta=true_log_mean_DR+delta)]
  sig=merge(rbind(tcl,tcr),sig,by=c("name","cat","id"))
  sig
}
ggplot(biases[name==name[1]&cat=="contact L"])+geom_point(aes(pos,etahat),alpha=0.1)+
  geom_line(aes(pos,eta))+facet_wrap(~ step)+xlim(550000,650000)+
  geom_line(aes(pos,true_eta),colour="red")





#bias around tad border
biases=cs@biases[,.(name,id,pos,true_log_rho,true_log_iota)]
biases=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  cbind(biases,cs@diagnostics$params[step==i&leg=="bias",.(step,log_rho=log_rho[[1]],log_iota=log_iota[[1]])])
}
ggplot(biases[name==name[1]])+geom_line(aes(pos,log_rho))+facet_wrap(~ step)+xlim(550000,650000)+
  geom_line(aes(pos,true_log_rho),colour="red")

#bias overall
biases=cs@biases[,.(name,id,pos,true_log_rho,true_log_iota)]
biases=foreach(i=1:cs@diagnostics$params[,max(step)],.combine=rbind) %do% {
  cbind(biases,cs@diagnostics$params[step==i&leg=="bias",.(step,log_rho=log_rho[[1]],log_iota=log_iota[[1]])])
}
ggplot(biases[name==name[1]])+geom_line(aes(pos,log_rho))+facet_wrap(~ step)+xlim(550000,650000)+
  geom_line(aes(pos,true_log_rho),colour="red")


#matrix zoom at tad border
signals[name==name[1]&unclass(bin1)>=29&unclass(bin1)<=32&unclass(bin2)>=29&unclass(bin2)<=32]
#signal
ggplot(cts[name==name[1]&unclass(bin)>=29&unclass(bin)<=32,.(phi=mean(phi)),by=c("bin","dbin")])+
  geom_point(aes(dbin,phi))+facet_wrap(~bin)
#mean with and w/0 signal
ggplot(cts[name==name[1]&unclass(bin)>=29&unclass(bin)<=32,.(lmu=mean(lmu),lsig=mean(lmu+phi)),by=c("bin","dbin")])+
  geom_point(aes(dbin,lmu,colour="w/o"))+geom_point(aes(dbin,lsig,colour="with"))+facet_wrap(~bin)
#weights in average z-score
a=cts[name==name[1]&unclass(bin)>=29&unclass(bin)<=32]
a[,c("z.bg","var.bg"):=list(count/exp(lmu)-1,(1/exp(lmu)+1/init$alpha))]
a[,c("z.sig","var.sig"):=list(count/exp(lmu+phi)-1,(1/exp(lmu+phi)+1/init$alpha))]
a=a[,.(weight.bg=mean(weight/var.bg),weight.sig=mean(weight/var.sig),
       zhat.bg=weighted.mean(z.bg, weight/var.bg), std.bg=1/sqrt(sum(weight/(2*var.bg))),
       zhat.sig=weighted.mean(z.sig, weight/var.sig), std.sig=1/sqrt(sum(weight/(2*var.sig)))),by=c("bin","dbin")]
ggplot(a)+geom_point(aes(dbin,weight.bg,colour="w/o"))+geom_point(aes(dbin,weight.sig,colour="with"))+facet_wrap(~bin)
#average z-score per dbin
ggplot(a)+geom_point(aes(dbin,zhat.bg,colour="w/o"))+geom_point(aes(dbin,zhat.sig,colour="with"))+facet_wrap(~bin)
#average std per dbin
ggplot(a)+geom_point(aes(dbin,std.bg,colour="w/o"))+geom_point(aes(dbin,std.sig,colour="with"))+facet_wrap(~bin)
#final z-scores for each bin
a=cts[unclass(bin)>=29&unclass(bin)<=32]
a[,c("z.bg","var.bg"):=list(count/exp(lmu)-1,(1/exp(lmu)+1/init$alpha))]
a[,c("z.sig","var.sig"):=list(count/exp(lmu+phi)-1,(1/exp(lmu+phi)+1/init$alpha))]
a=a[,.(weight.bg=mean(weight/var.bg),weight.sig=mean(weight/var.sig),
       zhat.bg=weighted.mean(z.bg, weight/var.bg), std.bg=1/sqrt(sum(weight/(2*var.bg))),
       zhat.sig=weighted.mean(z.sig, weight/var.sig), std.sig=1/sqrt(sum(weight/(2*var.sig)))),by=c("name","bin")]
ggplot(a)+geom_point(aes(bin,zhat.bg,colour="w/o"))+geom_point(aes(bin,zhat.sig,colour="with"))+facet_wrap(~ name)





a=Matrix(Diagonal(100))
a[60:70,60:70]=0.5
image(as.matrix(a))
b=Matrix(Diagonal(100))
b[60:70,]=-0.1
b[,60:70]=-0.1
b[60:70,60:70]=0.5
image(as.matrix(b))

#eigenvalues
ggplot(melt(data.table(x=seq(1,100),a=eigen(a)$values,b=eigen(b)$values),id.vars="x"))+geom_point(aes(x,value,colour=variable))
#eigenvectors
val=1
image(tcrossprod(eigen(a)$vectors[,val]))
image(tcrossprod(eigen(b,symmetric=T)$vectors[,val]))
ggplot(melt(data.table(x=seq(1,100),a=eigen(a)$vectors[,val], b=eigen(b,symmetric=T)$vectors[,val]), id.vars="x"))+geom_point(aes(x,value,colour=variable))



resolution=20000
ncores=30
cs=bin_all_datasets(cs, resolution=resolution, verbose=T, ncores=ncores)
mat=get_matrices(cs, resolution=resolution, group="all")
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=observed))+facet_wrap(~name)+
  scale_fill_gradient(high="red", low="white", na.value = "white")
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=expected))+facet_wrap(~name)+
  scale_fill_gradient(high="red", low="white", na.value = "white")


cs=detect_binned_interactions(cs, resolution=resolution, group="all", threshold=0.95, ncores=ncores)
save(cs,file="data/fake_csnorm_optimized.RData")
mat=get_interactions(cs, type="interactions", resolution=resolution, group="all", threshold=0.95, ref="expected")
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=signal))+facet_wrap(~name)+
  scale_fill_gradient(high="black", low="white", na.value = "white")+
  geom_point(aes(bin1,bin2,colour=direction),data=mat[is.significant==T])
mat[,.N,by=is.significant]

cs=group_datasets(cs, resolution=resolution, group="condition", ncores=ncores)
cs=detect_binned_interactions(cs, resolution=resolution, group="condition", threshold=0.95, ncores=ncores)
mat=get_interactions(cs, type="interactions", resolution=resolution, group="condition", threshold=0.95, ref="expected")
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=signal))+facet_wrap(~name)+
  scale_fill_gradient(high="black", low="white", na.value = "white")+
  geom_point(aes(bin1,bin2,colour=direction),data=mat[is.significant==T])
mat[,.N,by=is.significant]


cs=detect_binned_differences(cs, resolution=resolution, group="all", threshold=0.95,
                             ncores=ncores, ref=as.character(cs@experiments[1,name]))
save(cs,file="data/fake_csnorm_optimized.RData")
mat=get_interactions(cs, type="differences", resolution=resolution, group="all", threshold=0.95,
                     ref=as.character(cs@experiments[1,name]))
ggplot(mat)+geom_raster(aes(bin1,bin2,fill=difference))+facet_wrap(~name)+
  scale_fill_gradient(high="black", low="white", na.value = "white")+
  geom_point(aes(bin1,bin2,colour=direction),data=mat[is.significant==T])


cs=detect_binless_interactions(cs, resolution=resolution, group="all", ncores=30, niter=5, fit.decay=T)
mat=get_interactions(cs, type="binteractions", resolution=resolution, group="all", threshold=-1, ref="expected")
plot_binless_matrix(mat)


save(cs,file="data/fake_signal_csnorm_optimized.RData")





cs=detect_binless_differences(cs, resolution=resolution, group="all", ncores=ncores,
                              ref=as.character(cs@experiments[1,name]), niter=2)
save(cs, file=fname)


