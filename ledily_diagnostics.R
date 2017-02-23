library(csnorm)
library(foreach)
library(doParallel)
library(ggplot2)
library(data.table)

setwd("/home/yannick/simulations/cs_norm")

bpk=30
type="perf"

load("data/rao_HiCall_SELP_150k_csnorm_optimized_gauss_bpk30_perf_bpd10_bpb10_nostan.RData")
load("data/rao_HiCall_FOXP1_1.3M_csnorm_optimized_gauss_bpk3_nofill_perf_nostan.RData")
load(paste0("data/ledily_T47D_E60P30_gt150M_csnorm_optimized_bpk",bpk,"_",type,"_nostan.RData"))

bs=cs@par$biases
bs.new[,ori:="bpk30"]
bs=rbind(bs,bs.new)
ggplot(bs[cat=="contact L"&name=="GM MboI 1"])+
  geom_line(aes(pos,eta))+facet_grid(ori~.)

biases=rbind(cs@par$biases[,.(name,cat,run="15",etahat,pos,eta)],
             cs2@par$biases[,.(name,cat,run="16",etahat,pos,eta)])
  
ggplot(biases[name=="T47D es 60 MboI 1"])+facet_grid(run~cat)+
  geom_point(aes(pos,etahat,colour=cat),alpha=0.1)+
  geom_line(aes(pos,eta))

ggplot(cs@par$decay)+geom_point(aes(distance,kappahat,colour=name))+
  geom_line(aes(distance,kappa))+scale_x_log10()+facet_wrap(~name)

decays = foreach(i=cs2@diagnostics$params[leg=="decay",step], .combine=rbind) %do% {
  a=cs2@diagnostics$params[leg=="decay"&step==i,decay][[1]]
  a[,step:=i]
  a
}

id1=1140
id2=1352 
cs@biases[id==id1]
cs@biases[id==id2]
data=csd@data[begin1>149935719&begin1<150031830&begin2>149935719&begin2<150031830]
plot_raw(data, b1=149920719, b2=149920719, e1=150045830, e2=150045830,rsites = F)
plot_raw(data, b1=149935719, b2=149935719, e1=150031830, e2=150031830)
plot_raw(data, b1=149943719, b2=150019830, e1=149947719, e2=150023830)
plot_raw(data, b1=149985774, b2=149985774, e1=149989774, e2=149989774)
plot_raw(data, b1=149944719, b2=150020830, e1=149946719, e2=150022830)
cs@biases[pos>149944719&pos<149946719]
cs@biases[pos>150020830&pos<150022830]


ggplot(decays)+geom_point(aes(distance,kappahat,colour=name))+
  geom_line(aes(distance,kappa))+facet_grid(name~step)+scale_x_log10()
decays[,guys:=kappahat>kappa+2*std&distance>5e4&distance<1e5]
ggplot(decays[distance>7e4&distance<8e4])+geom_point(aes(distance,kappahat,colour=guys))+
  geom_line(aes(distance,kappa))+facet_grid(name~step)+scale_x_log10()
ggplot(decays)+geom_point(aes(distance,std,colour=name))+
  facet_grid(name~step)+scale_x_log10()+scale_y_log10()

decays[,rank:=abs(kappahat-kappa)/std]
decays[name=="T47D es 60 MboI 1"&step==10&distance>7e4&distance<8e4&guys==T][order(rank)]
decays[name=="T47D es 60 MboI 1"&step==10&distance>7e4&distance<8e4&guys==F][order(rank)]

dbin.good=decays[name=="T47D es 60 MboI 1"&step==10&distance>7e4&distance<8e4&guys==F][order(rank)][1,dbin]
dbin.bad=decays[name=="T47D es 60 MboI 1"&step==10&distance>7e4&distance<8e4&guys==T][order(rank)][.N,dbin]

biases = foreach(i=cs@diagnostics$params[leg=="bias",step], .combine=rbind) %do% {
  data.table(step=i,
             pos=1:length(cs@diagnostics$params[leg=="bias"&step==i,log_iota][[1]]),
             log_iota=cs@diagnostics$params[leg=="bias"&step==i,log_iota][[1]],
             log_rho=cs@diagnostics$params[leg=="bias"&step==i,log_rho][[1]])
}

ggplot(biases)+#[pos<1500&pos>1000])+
  geom_line(aes(pos,log_iota))+facet_grid(step~.)

ggplot(biases[,max(log_iota)-min(log_iota),by=step])+geom_line(aes(step,V1))

ggplot(counts[begin2<15100])+geom_raster(aes(begin1,begin2,fill=begin2-begin1>7.94328234724&begin2-begin1<8.12830516164))+facet_wrap(~name)
ggplot(counts[begin2<15100])+geom_raster(aes(begin1,begin2,fill=log(V1)))+facet_wrap(~name)
decays


zbias = csnorm:::get_nzeros_per_cutsite(cs, ncores=30)
csb = csnorm:::csnorm_gauss_genomic_muhat_mean(cs, zbias)
csb2 = csnorm:::csnorm_gauss_genomic_muhat_mean(cs2, zbias)
csb=rbind(csb$cts[,.(name,id,pos,cat,etahat,std,run="15")],
          csb2$cts[,.(name,id,pos,cat,etahat,std,run="16")])
setkey(csb,name,run,pos,cat)
ggplot(csb[name=="T47D es 60 MboI 1"])+geom_point(aes(pos,etahat,colour=cat),alpha=0.1)+facet_grid(cat~run)
ggplot(cts[name=="T47D es 60 MboI 1"&cat=="contact L"])+
  geom_point(aes(pos,eta),alpha=0.1)+facet_wrap(~weight>1)
