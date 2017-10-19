library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

resolution=5000
sub="Peak1_450k"

load(paste0("data/rao_HiCall_GM12878_",sub,"_csnorm_optimized_base5k_bpk",30,"_qmin",0.0,"_bpb",10,"_bpd",10,".RData"))


#matrix plots
mat=get_interactions(cs, type="binteractions", resolution=resolution, group="all", threshold=-1, ref="expected")
ggplot(mat)+geom_raster(aes(begin2,begin1,fill=pmin(phi,3)))+
  geom_raster(aes(begin1,begin2,fill=pmin(phihat,3)))+
  facet_wrap(~name)+scale_fill_gradient(high=muted("red"), low="white", na.value = "white")+stat_function(fun=identity)+
  geom_point(data=cs@biases[pos-shift(pos)==175,.(pos)],aes(pos,pos))

ggplot(mat[begin1==begin2])+geom_point(aes(begin1,phi,colour=name))

ggplot(mat[begin2<74000000])+geom_raster(aes(begin1,begin2,fill=phihat))+
  geom_raster(aes(begin2,begin1,fill=phi))+
  facet_wrap(~name)+scale_fill_gradient(high=muted("red"), low="white", na.value = "white")+stat_function(fun=identity)+
  stat_function(fun=function(x){73885158})

ggplot(mat[bin1==bin2])+geom_pointrange(aes(bin1,phihat,ymin=phihat-sqrt(phihat.var),ymax=phihat+sqrt(phihat.var)))

ggplot(mat[bin1==bin2])+geom_point(aes(bin1,weight,colour=phihat))+scale_colour_gradient2()
ggplot(mat[bin1==bin2])+geom_point(aes(ncounts,weight,colour=phihat))+scale_colour_gradient2()


#count sums
bts=melt(cs@biases[name==name[1]],id.vars=c("name","id","pos"))[,.(id,pos,variable,value)]
cts=rbind(cs@counts[name==name[1],.(id=id1, pos=pos1, contact.R=contact.close, contact.L=contact.far)],
          cs@counts[name==name[1],.(id=id1, pos=pos1, contact.R=contact.down,  contact.L=contact.up)],
          cs@counts[name==name[1],.(id=id2, pos=pos2, contact.R=contact.far,   contact.L=contact.close)],
          cs@counts[name==name[1],.(id=id2, pos=pos2, contact.R=contact.down,  contact.L=contact.up)])
cts=rbind(bts,melt(cts[,.(contact.R=sum(contact.R),contact.L=sum(contact.L)),
                       by=c("id","pos")],id.vars=c("id","pos")))
cts=cts[,.(value=sum(value)),by="pos"]

ggplot(cts[,.(dist=pos-shift(pos),value)][dist<500])+geom_point(aes(dist,value))#+scale_x_log10()
ggplot(cts[,.(dist=pos-shift(pos),value)][dist<250&dist>100&dist!=175])+geom_point(aes(dist,value))#+scale_x_log10()
cts[,.(dist=pos-shift(pos),value)][dist==175]

ggplot(cs@biases[name==name[1],.(dist=pos-shift(pos))][dist<500])+geom_histogram(aes(dist),bins=100)#+scale_x_log10()
ggplot(b[,.(dist=pos-shift(pos),ori)][dist<500])+geom_histogram(aes(dist),bins=100)+facet_grid(ori~.)


ggplot(cs@biases[,.(pos,var=pos-shift(pos)==175)])+geom_point(aes(pos,1+var))

#vary resolutions
resolution=5000
#zeros and cts
zeros = binless:::get_nzeros_binning(cs, resolution, ncores = 30)
cts = binless:::predict_binned_counts_irls(cs, resolution, zeros)
setkeyv(cts,c("name","bin1","bin2"))
#phi and trails
names=cs@experiments[,unique(name)]
groups=data.table(name=names,groupname=names)
setkey(groups,name)
stuff = binless:::prepare_signal_matrix(cs, groups, resolution)
mat=stuff$mat
trails=stuff$trails
#phihat
mat = binless:::compute_raw_signal(cts, cs@par$alpha, mat)
bin1.begin=mat[,bin1]
bin2.begin=mat[,bin2]
levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
mat[,begin1:=as.integer(as.character(bin1.begin))]
mat[,begin2:=as.integer(as.character(bin2.begin))]
mat[,res:=resolution]

#plots
ggplot(mat[begin2<73910000&begin1>73870000])+geom_raster(aes(begin1,begin2,fill=pmin(phihat,3)))+
  geom_raster(aes(begin2,begin1,fill=phi))+
  facet_wrap(~name)+scale_fill_gradient(high=muted("red"), low="white", na.value = "white")+stat_function(fun=identity)+
  stat_function(fun=function(x){73888158})

ggplot(rbind(mat,ref.mat)[bin1==bin2&begin2<74000000])+geom_pointrange(aes(begin1,phihat,colour=factor(res),
                                                           ymin=phihat-sqrt(phihat.var),ymax=phihat+sqrt(phihat.var)))
ggplot(rbind(mat,ref.mat)[bin1==bin2&begin2<73900000&begin1>73860000])+geom_pointrange(aes(begin1,phihat,colour=factor(res),
                                                                           ymin=phihat-sqrt(phihat.var),ymax=phihat+sqrt(phihat.var)))

ggplot(ref.mat[bin1==bin2][bin1==bin2&begin2<73900000&begin1>73860000])+geom_point(aes(begin1,weight,colour=phihat))+scale_colour_gradient2()
ggplot(rbind(mat,ref.mat)[bin1==bin2][bin1==bin2&begin2<73900000&begin1>73860000])+geom_point(aes(begin1,weight,colour=factor(res)))#+scale_colour_gradient2()
ggplot(mat[unclass(bin1)==unclass(bin2)-1&begin2<73890158&begin1>73887158])+geom_point(aes(begin1,weight,colour=factor(res)))#+scale_colour_gradient2()

ggplot(ref.mat[bin1==bin2][bin1==bin2&begin2<73900000&begin1>73860000])+geom_point(aes(ncounts,weight,colour=phihat))+scale_colour_gradient2()
ggplot(mat[bin1==bin2][bin1==bin2&begin2<73900000&begin1>73860000])+geom_point(aes(ncounts,weight,colour=phihat))+scale_colour_gradient2()
ggplot(rbind(mat,ref.mat)[bin1==bin2][bin1==bin2&begin2<73900000&begin1>73860000])+geom_point(aes(ncounts,weight,colour=factor(res)))#+scale_colour_gradient2()

bin1.begin=zeros[,bin1]
bin2.begin=zeros[,bin2]
levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
zeros[,begin1:=as.integer(as.character(bin1.begin))]
zeros[,begin2:=as.integer(as.character(bin2.begin))]
zeros[,res:=resolution]

ggplot(zeros[begin2<73910000&begin1>73870000])+geom_raster(aes(begin1,begin2,fill=nnz))+
  geom_raster(aes(begin2,begin1,fill=nnz))+
  facet_wrap(~name)+scale_fill_gradient(high=muted("red"), low="white", na.value = "white")+stat_function(fun=identity)+
  stat_function(fun=function(x){73888158})


sbins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
ncuts=cs@biases[,.(name,pos,bin=cut(pos, sbins,
                       ordered_result=T, right=F, include.lowest=T,dig.lab=12))][,.N,by=c("name","bin")]
bin1.begin=ncuts[,bin]
levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
ncuts[,begin:=as.integer(as.character(bin1.begin))]
ncuts[,res:=resolution]

ggplot(ncuts[begin<73910000&begin>73870000])+geom_point(aes(begin,N,colour=begin==73888158))
ggplot(ncuts[begin<73890000&begin>73880000])+geom_point(aes(begin,N,colour=begin==73885158))

ggplot(cs@biases[pos<73890000&pos>73880000,.(name,pos,bin=cut(pos, sbins, ordered_result=T,
                                     right=F, include.lowest=T,dig.lab=12))])+
  geom_density(aes(pos))+
  geom_bar(aes(pos,5e-5,colour=bin),stat="identity")+guides(colour=F)


ggplot(cs@biases[pos<73890000&pos>73880000,.(name,pos,bin=cut(pos, sbins,
                              ordered_result=T, right=F, include.lowest=T,dig.lab=12))][
                                ,mean(pos-shift(pos),na.rm=T),by=c("name","bin")])+
  geom_point(aes(bin,V1))

ggplot(cs@biases[pos<73890000&pos>73880000,.(name,dist=pos-shift(pos),bin=cut(pos, sbins,
                                                              ordered_result=T, right=F, include.lowest=T,dig.lab=12))])+
  geom_boxplot(aes(bin,dist))

ggplot(cts[pos<73890000&pos>73880000,.(pos,value,bin=cut(pos, sbins,
                                                  ordered_result=T, right=F, include.lowest=T,dig.lab=12))])+
  geom_jitter(aes(bin,value,colour=bin))

cts[pos<73890000&pos>73880000,
    .(pos,value,bin=cut(pos, sbins,ordered_result=T, right=F, include.lowest=T,dig.lab=12))][
      ,sum(value),by=bin]


ggplot(cs@biases[,.(name,pos,dist=pos-shift(pos))])+geom_histogram(aes(dist),bins=100)+scale_x_log10()
ggplot(cs@biases[,.(name,pos,dist=pos-shift(pos))][dist<1000])+geom_histogram(aes(dist),bins=100)
ggplot(cs@biases[,.(name,pos,dist=pos-shift(pos))][dist<100])+geom_histogram(aes(dist),bins=100)

load("data/rao_HiCall_IMR90_SEMA3C_1M_csdata.RData")
load("data/rao_HiCall_GM12878_SEMA3C_1M_csdata.RData")
load("data/ledily_T47D_es_60_csdata.RData")
load("data/ledily_T47D_pg_180_csdata.RData")
load("data/ralph_EScell_Rbfox1_3.5M_csdata.RData")
ggplot(csd@biases[,.(dist=pos-shift(pos))][dist<500])+geom_histogram(aes(dist),bins=100)#+scale_x_log10()

bts=melt(csd@biases,id.vars=c("name","id","pos"))[,.(id,pos,variable,value)]
cts=rbind(csd@counts[,.(id=id1, pos=pos1, contact.R=contact.close, contact.L=contact.far)],
          csd@counts[,.(id=id1, pos=pos1, contact.R=contact.down,  contact.L=contact.up)],
          csd@counts[,.(id=id2, pos=pos2, contact.R=contact.far,   contact.L=contact.close)],
          csd@counts[,.(id=id2, pos=pos2, contact.R=contact.down,  contact.L=contact.up)])
cts=rbind(bts,melt(cts[,.(contact.R=sum(contact.R),contact.L=sum(contact.L)),
                       by=c("id","pos")],id.vars=c("id","pos")))
cts=cts[,.(value=sum(value)),keyby="pos"]
ggplot(cts[,.(dist=pos-shift(pos),value)][dist<500])+geom_point(aes(dist,value))+geom_smooth(aes(dist,value))
ggplot(cts[,.(dist=pos-shift(pos),value)][dist<500&value<20])+geom_point(aes(dist,value))#+scale_x_log10()




load("data/caulo_BglIIr1_all_csdata.RData")
load("data/caulo_NcoI_all_csdata.RData")
load("data/rao_HiC035_HindIII_GM12878_Peak1_450k_csdata.RData")
load("data/rao_HiC036_NcoI_GM12878_Peak1_450k_csdata.RData")
ggplot(csd@biases[,.(dist=pos-shift(pos))][dist<5000])+geom_histogram(aes(dist),bins=100)#+scale_x_log10()



mat=mat[,.(name,bin1,bin2,phi=0)]
cts.cp = mat[cts,,on=c("name","bin1","bin2")]
cts.cp[,c("z","var"):=list(count/exp(phi+eC+lmu.base+log_decay)-1,
                           (1/exp(phi+eC+lmu.base+log_decay)+1/dispersion))]
mat = cts.cp[,.(phihat=weighted.mean(z+phi, weight/var),
                phihat.var=1/sum(weight/var),
                ncounts=sum(weight)),keyby=c("name","bin1","bin2")][mat,,on=c("name","bin1","bin2")]
mat[is.na(phihat),c("phihat","phihat.var","ncounts"):=list(1,Inf,0)] #bins with no detectable counts
mat[,c("valuehat","weight"):=list(phihat,1/phihat.var)]
bin1.begin=mat[,bin1]
bin2.begin=mat[,bin2]
levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
mat[,begin1:=as.integer(as.character(bin1.begin))]
mat[,begin2:=as.integer(as.character(bin2.begin))]
mat[,res:=resolution]
ggplot(mat[bin1==bin2&begin2<73900000&begin1>73860000])+geom_pointrange(aes(begin1,phihat,colour=factor(res),
                                                                                           ymin=phihat-sqrt(phihat.var),ymax=phihat+sqrt(phihat.var)))

bin1.begin=cts.cp[,bin1]
bin2.begin=cts.cp[,bin2]
levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
cts.cp[,begin1:=as.integer(as.character(bin1.begin))]
cts.cp[,begin2:=as.integer(as.character(bin2.begin))]
cts.cp[,res:=resolution]
#weight
ggplot(cts.cp[bin1==bin2&begin2<73900000&begin1>73860000])+geom_jitter(aes(begin1,weight,colour=factor(count)))
#var
ggplot(cts.cp[bin1==bin2&begin2<73900000&begin1>73860000])+geom_jitter(aes(begin1,var,colour=factor(count)))
#weight/var
ggplot(cts.cp[bin1==bin2&begin2<73900000&begin1>73860000])+geom_jitter(aes(begin1,weight/var,colour=factor(count)))
#z
cts.cp[,alpha:=weight/var]
cts.cp[,alpha:=alpha/max(alpha),by=c("bin1","bin2")]
ggplot(cts.cp[bin1==bin2&begin2<73900000&begin1>73860000])+geom_jitter(aes(begin1,z,colour=factor(count),alpha=alpha))
