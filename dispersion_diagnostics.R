library(binless)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)
library(methods)
library(igraph)
library(MASS)

if (F) {
  #load a normalized dataset, then do this
  counts=binless:::fill_zeros(cs@counts,cs@biases,dmin=cs@settings$dmin)
  counts=binless:::predict_all_parallel(cs,counts,ncores = ncores)
  setkeyv(counts,key(cs@counts))
  stopifnot(cs@biases[,.N]==length(cs@par$log_iota))
  #add signal contribution if available
  if (length(cs@settings$sbins)>2) {
    counts[,bin1:=cut(pos1, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    counts[,bin2:=cut(pos2, cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    counts=cs@par$signal[,.(name,bin1,bin2,phi)][counts,,on=key(cs@par$signal)]
    counts[,c("log_mean_cclose","log_mean_cfar","log_mean_cup","log_mean_cdown"):=
             list(log_mean_cclose+phi,log_mean_cfar+phi,log_mean_cup+phi,log_mean_cdown+phi)]
  }
  biases=cs@biases
  biases[,log_iota:=cs@par$log_iota]
  biases[,log_rho:=cs@par$log_rho]
  biases[,log_mean_DL:=cs@par$eDE[1]+log_iota]
  biases[,log_mean_DR:=cs@par$eDE[1]+log_rho]
  biases[,log_mean_RJ:=cs@par$eRJ[1]+(log_rho+log_iota)/2.]
}

#load("tmp/dispersion_diag.RData")

counts=rbind(counts[,.(name,bin1,bin2,phi,id1,id2,pos1,pos2,distance,log_decay,count=contact.close,log_mean=log_mean_cclose,type="contact",subtype="close")],
             counts[,.(name,bin1,bin2,phi,id1,id2,pos1,pos2,distance,log_decay,count=contact.down,log_mean=log_mean_cdown,type="contact",subtype="down")],
             counts[,.(name,bin1,bin2,phi,id1,id2,pos1,pos2,distance,log_decay,count=contact.far,log_mean=log_mean_cfar,type="contact",subtype="far")],
             counts[,.(name,bin1,bin2,phi,id1,id2,pos1,pos2,distance,log_decay,count=contact.up,log_mean=log_mean_cup,type="contact",subtype="up")])
counts[,c("type","subtype"):=list(factor(type),factor(subtype))]
biases=rbind(biases[,.(name,id,pos,count=dangling.L,log_mean=log_mean_DL,type="dangling",subtype="DL")],
             biases[,.(name,id,pos,count=dangling.R,log_mean=log_mean_DR,type="dangling",subtype="DR")],
             biases[,.(name,id,pos,count=rejoined,log_mean=log_mean_RJ,type="rejoined",subtype="RJ")])
biases[,c("type","subtype"):=list(factor(type),factor(subtype))]



### basic model: fit all counts with glm.nb: theta is 0.7715 +- 0.0137
alldata=rbind(counts[,.(name,count,log_mean,type,subtype)],biases[,.(name,count,log_mean,type,subtype)])
fit=glm.nb(count~log_mean+type,link=log,init.theta=cs@par$alpha,data=alldata)
summary(fit)


### basic model: same with stan: theta is 0.76659
data = list( ncounts=alldata[,.N], count=alldata[,count], log_mean=alldata[,log_mean],
             ntypes=alldata[,nlevels(type)], type=alldata[,unclass(type)],
             ngroups=1, dispgroup=alldata[,rep(1,.N)])
stanmodel = stan_model("dispersion_diagnostics.stan")
op1=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)


### one dispersion for each count type
data = list( ncounts=alldata[,.N], count=alldata[,count], log_mean=alldata[,log_mean],
             ntypes=alldata[,nlevels(type)], type=alldata[,unclass(type)],
             ngroups=alldata[,nlevels(type)], dispgroup=alldata[,unclass(type)])
stanmodel = stan_model("dispersion_diagnostics.stan")
op2=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
data.table(name=alldata[,levels(type)], exposure=op2$par$exposure, dispersion=op2$par$dispersion)
"
       name     exposure dispersion
1:  contact  0.068471138  0.6608148
2: dangling -0.001763415  0.8802533
3: rejoined -0.014559066  1.8779235
"


### one dispersion for each count subtype
data = list( ncounts=alldata[,.N], count=alldata[,count], log_mean=alldata[,log_mean],
             ntypes=alldata[,nlevels(type)], type=alldata[,unclass(type)],
             ngroups=alldata[,nlevels(subtype)], dispgroup=alldata[,unclass(subtype)])
stanmodel = stan_model("dispersion_diagnostics.stan")
op3=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
data.table(name=alldata[,levels(type)], exposure=op3$par$exposure)
"
name      exposure
1:  contact  0.0684446953
2: dangling -0.0003403972
3: rejoined -0.0147130569
"
data.table(name=alldata[,levels(subtype)], dispersion=op3$par$dispersion)
'    name dispersion
1: close  0.5966946
2:  down  0.8758821
3:   far  0.6605046
4:    up  0.5568410
5:    DL  0.8459249
6:    DR  0.9176570
7:    RJ  1.8712655
'


### group by log10-distance, one exposure for all
ngroups=10
dbins=seq(counts[,min(log10(distance))],counts[,max(log10(distance))],length.out = ngroups+1)
counts[,dbin:=cut(log10(distance),dbins, include.lowest=T, ordered_result=T,dig.lab=12,right=F)]
counts[,.N,by=dbin][,min(N)]
#plot decay
ggplot(cs@par$decay)+geom_line(aes(distance,kappa))+geom_point(aes(distance,kappahat))+scale_x_log10()
#ggplot(counts[,.(distance=10^mean(log10(distance)),count=mean(count),std=sd(count),mu=mean(exp(log_mean))),by=dbin])+
#  geom_line(aes(distance,mu))+geom_point(aes(distance,count))+scale_x_log10()
data = list( ncounts=counts[,.N], count=counts[,count], log_mean=counts[,log_mean],
             ntypes=counts[,nlevels(type)], type=counts[,unclass(type)],
             ngroups=counts[,nlevels(dbin)], dispgroup=counts[,unclass(dbin)])
op4=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
data.table(name=counts[,levels(type)], exposure=op4$par$exposure)
"
      name   exposure
1: contact 0.06857623
"
data.table(name=counts[,levels(dbin)], distance=counts[,10^mean(log10(distance)),by=dbin]$V1, dispersion=op4$par$dispersion)
"
                             name   distance dispersion
 1:             [3,3.28449702145)   1439.828 1.61808452
 2: [3.28449702145,3.56899404291)   2756.920 1.21073636
 3: [3.56899404291,3.85349106436)   5337.913 1.15518920
 4: [3.85349106436,4.13798808582)  10267.700 1.08679261
 5: [4.13798808582,4.42248510727)  19731.051 0.69248180
 6: [4.42248510727,4.70698212873)  38026.981 0.44236278
 7: [4.70698212873,4.99147915018)  72947.305 0.25827075
 8: [4.99147915018,5.27597617163) 139626.572 0.17159361
 9: [5.27597617163,5.56047319309) 265689.873 0.11475870
10: [5.56047319309,5.84497021454] 469219.142 0.09917161
"
ggplot(data.table(distance=counts[,10^mean(log10(distance)),by=dbin]$V1, dispersion=op4$par$dispersion))+geom_point(aes(distance,dispersion))+scale_x_log10()
ggplot(data.table(count=counts[,mean(count),by=dbin]$V1, dispersion=op4$par$dispersion))+geom_point(aes(count,dispersion))+scale_x_log10()



### group by log10-distance, one exposure for every group
data = list( ncounts=counts[,.N], count=counts[,count], log_mean=counts[,log_mean],
             ntypes=counts[,nlevels(dbin)], type=counts[,unclass(dbin)],
             ngroups=counts[,nlevels(dbin)], dispgroup=counts[,unclass(dbin)])
op5=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
data.table(name=counts[,levels(dbin)], distance=counts[,10^mean(log10(distance)),by=dbin]$V1, exposure=op5$par$exposure, dispersion=op5$par$dispersion)
"
name   distance   exposure dispersion
1:             [3,3.28449702145)   1439.828 0.10348742 1.63233539
2: [3.28449702145,3.56899404291)   2756.920 0.10380416 1.22234060
3: [3.56899404291,3.85349106436)   5337.913 0.09026299 1.16348002
4: [3.85349106436,4.13798808582)  10267.700 0.08177860 1.09094157
5: [4.13798808582,4.42248510727)  19731.051 0.08553965 0.69452808
6: [4.42248510727,4.70698212873)  38026.981 0.08204920 0.44351309
7: [4.70698212873,4.99147915018)  72947.305 0.05695900 0.25785053
8: [4.99147915018,5.27597617163) 139626.572 0.04027585 0.17114312
9: [5.27597617163,5.56047319309) 265689.873 0.04704953 0.11477032
10: [5.56047319309,5.84497021454] 469219.142 0.03824811 0.09861607
"
ggplot(data.table(distance=counts[,10^mean(log10(distance)),by=dbin]$V1, d4=op4$par$dispersion, d5=op5$par$dispersion))+
  geom_point(aes(distance,d5),colour="green")+geom_point(aes(distance,d4),colour="red")+scale_x_log10()


### group by linear distance, one exposure for every group
ngroups=10
dbins=seq(counts[,min((distance))],counts[,max((distance))],length.out = ngroups+1)
counts[,dbin:=cut((distance),dbins, include.lowest=T, ordered_result=T,dig.lab=12,right=F)]
counts[,.N,by=dbin][,min(N)]
data = list( ncounts=counts[,.N], count=counts[,count], log_mean=counts[,log_mean],
             ntypes=counts[,nlevels(dbin)], type=counts[,unclass(dbin)],
             ngroups=counts[,nlevels(dbin)], dispgroup=counts[,unclass(dbin)])
op6=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
data.table(name=counts[,levels(dbin)], distance=counts[,10^mean(log10(distance)),by=dbin]$V1, exposure=op6$par$exposure, dispersion=op6$par$dispersion)
"
                   name  distance      exposure dispersion
 1:      [1000,70879.4)  27095.31  0.0866708717 0.93382956
 2:  [70879.4,140758.8) 103199.21  0.0378065814 0.20224480
 3: [140758.8,210638.2) 173687.24  0.0418016315 0.17018825
 4: [210638.2,280517.6) 243897.17  0.0391916209 0.11750810
 5:   [280517.6,350397) 313775.46  0.0624229141 0.09441897
 6:   [350397,420276.4) 383567.27  0.0612794884 0.08846093
 7: [420276.4,490155.8) 453271.74 -0.0006303332 0.08190514
 8: [490155.8,560035.2) 522295.23  0.0821628176 0.17331689
 9: [560035.2,629914.6) 590869.04  0.0961419703 0.07392706
10:   [629914.6,699794] 652803.16 -0.5628318701 0.02511213
"
ggplot(data.table(distance=counts[,mean(as.numeric(distance)),by=dbin]$V1, dispersion=op6$par$dispersion))+geom_point(aes(distance,dispersion))
ggplot(data.table(count=counts[,mean(count),by=dbin]$V1, dispersion=op6$par$dispersion))+geom_point(aes(count,dispersion))+scale_x_log10()





### group by ligation efficiency
gdata=rbind(counts[,.(id=id1,count)],counts[,.(id=id2,count)])[,.(count=sum(count)/2),by=id]
ggplot(gdata)+geom_histogram(aes(count))#+scale_x_log10()
gbins=quantile(gdata[,count],c(0,1/3,2/3,1))
gdata[,gbin:=cut2(count,g=3)]
gcounts=merge(counts,gdata[,.(id1=id,gbin1=gbin)],by="id1")
gcounts=merge(gcounts,gdata[,.(id2=id,gbin2=gbin)],by="id2")
gcounts[,gbin:=factor(unclass(gbin1)*unclass(gbin2))]
data = list( ncounts=gcounts[,.N], count=gcounts[,count], log_mean=gcounts[,log_mean],
             ntypes=gcounts[,nlevels(gbin)], type=gcounts[,unclass(gbin)],
             ngroups=gcounts[,nlevels(gbin)], dispgroup=gcounts[,unclass(gbin)])
op7=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
data.table(name=gcounts[,levels(gbin)], count=gcounts[,mean(count),by=gbin]$V1, exposure=op7$par$exposure, dispersion=op7$par$dispersion)
"
   name      count   exposure dispersion
1:    1 0.00334087 -0.8668380  0.1770302
2:    2 0.00677934 -0.3609118  0.3215196
3:    3 0.01474459 -0.2216776  0.3971618
4:    4 0.01066113  0.1665895  0.5118275
5:    6 0.02417419  0.2736153  0.7247421
6:    9 0.04310824  0.3593629  1.1305923
"
ggplot(data.table(count=gcounts[,mean(count),by=gbin]$V1, dispersion=op7$par$dispersion))+geom_point(aes(count,dispersion))
ggplot(gcounts[,.N,by=c("bin1","bin2","gbin")][,.SD[N==max(N)],keyby=c("bin1","bin2")])+geom_raster(aes(bin1,bin2,fill=gbin))


#group by ligation efficiency AND log10-distance bin
ggplot(gcounts[,.N,by=c("bin1","bin2","gbin")][,.SD[N==max(N)],keyby=c("bin1","bin2")])+geom_raster(aes(bin1,bin2,fill=gbin))
ggplot(gcounts[,.N,by=c("bin1","bin2","dbin")][,.SD[N==max(N)],keyby=c("bin1","bin2")])+geom_raster(aes(bin1,bin2,fill=dbin))
gcounts[,mixbin:=interaction(gbin,dbin)]
ggplot(gcounts[,.N,by=c("bin1","bin2","mixbin")][,.SD[N==max(N)],keyby=c("bin1","bin2")])+geom_raster(aes(bin1,bin2,fill=mixbin))
data = list( ncounts=gcounts[,.N], count=gcounts[,count], log_mean=gcounts[,log_mean],
             ntypes=gcounts[,nlevels(mixbin)], type=gcounts[,unclass(mixbin)],
             ngroups=gcounts[,nlevels(mixbin)], dispgroup=gcounts[,unclass(mixbin)])
op8=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
result=cbind(gcounts[,.(count=mean(count),mu=mean(exp(log_mean))),keyby=c("mixbin","gbin","dbin")],exposure=op8$par$exposure,dispersion=op8$par$dispersion)
"
                             mixbin gbin                          dbin       count    exposure dispersion
 1:             1.[3,3.28449702145)    1             [3,3.28449702145) 0.043650794 -0.95384379 5.75516999
2:             2.[3,3.28449702145)    2             [3,3.28449702145) 0.108946213 -0.28022195 1.47818969
3:             3.[3,3.28449702145)    3             [3,3.28449702145) 0.198407202 -0.05610748 1.10267898
4:             4.[3,3.28449702145)    4             [3,3.28449702145) 0.213400901  0.17132181 1.51507705
5:             6.[3,3.28449702145)    6             [3,3.28449702145) 0.352380952  0.29051863 1.98741310
6:             9.[3,3.28449702145)    9             [3,3.28449702145) 0.567256637  0.34687581 2.16681498
7: 1.[3.28449702145,3.56899404291)    1 [3.28449702145,3.56899404291) 0.048087141 -0.77955637 0.31934747
8: 2.[3.28449702145,3.56899404291)    2 [3.28449702145,3.56899404291) 0.092419825 -0.32538685 0.98702742
9: 3.[3.28449702145,3.56899404291)    3 [3.28449702145,3.56899404291) 0.157053509 -0.17834719 0.67062627
10: 4.[3.28449702145,3.56899404291)    4 [3.28449702145,3.56899404291) 0.198019802  0.23310337 0.82527698
11: 6.[3.28449702145,3.56899404291)    6 [3.28449702145,3.56899404291) 0.320312500  0.28188137 1.42031881
12: 9.[3.28449702145,3.56899404291)    9 [3.28449702145,3.56899404291) 0.549056604  0.38438498 2.04987963
13: 1.[3.56899404291,3.85349106436)    1 [3.56899404291,3.85349106436) 0.036299127 -0.68128566 0.33719324
14: 2.[3.56899404291,3.85349106436)    2 [3.56899404291,3.85349106436) 0.066939469 -0.31692137 0.60962065
15: 3.[3.56899404291,3.85349106436)    3 [3.56899404291,3.85349106436) 0.106182326 -0.24851726 0.80838038
16: 4.[3.56899404291,3.85349106436)    4 [3.56899404291,3.85349106436) 0.143378063  0.20500554 0.94883061
17: 6.[3.56899404291,3.85349106436)    6 [3.56899404291,3.85349106436) 0.235443221  0.27719558 1.50445205
18: 9.[3.56899404291,3.85349106436)    9 [3.56899404291,3.85349106436) 0.397445255  0.34286642 1.59766807
19: 1.[3.85349106436,4.13798808582)    1 [3.85349106436,4.13798808582) 0.014659271 -1.01292232 0.35099103
20: 2.[3.85349106436,4.13798808582)    2 [3.85349106436,4.13798808582) 0.038746259 -0.31395569 0.37210194
21: 3.[3.85349106436,4.13798808582)    3 [3.85349106436,4.13798808582) 0.064266901 -0.21923177 0.59846816
22: 4.[3.85349106436,4.13798808582)    4 [3.85349106436,4.13798808582) 0.079611949  0.13923013 0.70892512
23: 6.[3.85349106436,4.13798808582)    6 [3.85349106436,4.13798808582) 0.141880473  0.29584129 1.21470550
24: 9.[3.85349106436,4.13798808582)    9 [3.85349106436,4.13798808582) 0.232268544  0.31387769 2.01973309
25: 1.[4.13798808582,4.42248510727)    1 [4.13798808582,4.42248510727) 0.010753486 -0.79696866 0.11928730
26: 2.[4.13798808582,4.42248510727)    2 [4.13798808582,4.42248510727) 0.021316072 -0.36826236 0.26207474
27: 3.[4.13798808582,4.42248510727)    3 [4.13798808582,4.42248510727) 0.039620170 -0.14805745 0.47292421
28: 4.[4.13798808582,4.42248510727)    4 [4.13798808582,4.42248510727) 0.045878631  0.15399849 0.66280488
29: 6.[4.13798808582,4.42248510727)    6 [4.13798808582,4.42248510727) 0.080117484  0.30507063 0.69798865
30: 9.[4.13798808582,4.42248510727)    9 [4.13798808582,4.42248510727) 0.124231154  0.31054687 1.17760609
31: 1.[4.42248510727,4.70698212873)    1 [4.42248510727,4.70698212873) 0.006128181 -0.78394202 0.10560625
32: 2.[4.42248510727,4.70698212873)    2 [4.42248510727,4.70698212873) 0.012273062 -0.32811025 0.17447008
33: 3.[4.42248510727,4.70698212873)    3 [4.42248510727,4.70698212873) 0.019809807 -0.22747205 0.25139629
34: 4.[4.42248510727,4.70698212873)    4 [4.42248510727,4.70698212873) 0.026533275  0.20197217 0.34131088
35: 6.[4.42248510727,4.70698212873)    6 [4.42248510727,4.70698212873) 0.041051681  0.25962178 0.40064227
36: 9.[4.42248510727,4.70698212873)    9 [4.42248510727,4.70698212873) 0.069740611  0.38799683 0.95687751
37: 1.[4.70698212873,4.99147915018)    1 [4.70698212873,4.99147915018) 0.002868187 -1.01615322 0.08616136
38: 2.[4.70698212873,4.99147915018)    2 [4.70698212873,4.99147915018) 0.006832209 -0.38607303 0.11203586
39: 3.[4.70698212873,4.99147915018)    3 [4.70698212873,4.99147915018) 0.010363216 -0.33887287 0.16666587
40: 4.[4.70698212873,4.99147915018)    4 [4.70698212873,4.99147915018) 0.014651537  0.14800136 0.15808669
41: 6.[4.70698212873,4.99147915018)    6 [4.70698212873,4.99147915018) 0.024512382  0.29251081 0.29970935
42: 9.[4.70698212873,4.99147915018)    9 [4.70698212873,4.99147915018) 0.040033844  0.38729586 0.40172958
43: 1.[4.99147915018,5.27597617163)    1 [4.99147915018,5.27597617163) 0.002289577 -0.82893675 0.03015193
44: 2.[4.99147915018,5.27597617163)    2 [4.99147915018,5.27597617163) 0.004354653 -0.40371023 0.05452904
45: 3.[4.99147915018,5.27597617163)    3 [4.99147915018,5.27597617163) 0.007214537 -0.25468188 0.10617965
46: 4.[4.99147915018,5.27597617163)    4 [4.99147915018,5.27597617163) 0.009318914  0.14186523 0.13345303
47: 6.[4.99147915018,5.27597617163)    6 [4.99147915018,5.27597617163) 0.014521628  0.23520573 0.16828700
48: 9.[4.99147915018,5.27597617163)    9 [4.99147915018,5.27597617163) 0.023832246  0.36105664 0.34454358
49: 1.[5.27597617163,5.56047319309)    1 [5.27597617163,5.56047319309) 0.001359536 -0.98172314 0.02123503
50: 2.[5.27597617163,5.56047319309)    2 [5.27597617163,5.56047319309) 0.002999643 -0.40917895 0.07417872
51: 3.[5.27597617163,5.56047319309)    3 [5.27597617163,5.56047319309) 0.005185686 -0.20363250 0.07150842
52: 4.[5.27597617163,5.56047319309)    4 [5.27597617163,5.56047319309) 0.006525350  0.15609167 0.10962722
53: 6.[5.27597617163,5.56047319309)    6 [5.27597617163,5.56047319309) 0.010225732  0.25751490 0.12679390
54: 9.[5.27597617163,5.56047319309)    9 [5.27597617163,5.56047319309) 0.016532454  0.39974095 0.17154545
55: 1.[5.56047319309,5.84497021454]    1 [5.56047319309,5.84497021454] 0.001027131 -0.85968017 0.01693903
56: 2.[5.56047319309,5.84497021454]    2 [5.56047319309,5.84497021454] 0.002222855 -0.35777483 0.04331523
57: 3.[5.56047319309,5.84497021454]    3 [5.56047319309,5.84497021454] 0.003834840 -0.18012006 0.07957250
58: 4.[5.56047319309,5.84497021454]    4 [5.56047319309,5.84497021454] 0.004982612  0.16564625 0.09445443
59: 6.[5.56047319309,5.84497021454]    6 [5.56047319309,5.84497021454] 0.008026975  0.25883524 0.13449538
60: 9.[5.56047319309,5.84497021454]    9 [5.56047319309,5.84497021454] 0.011920371  0.32652878 0.11519065
"
ggplot(result)+geom_line(aes(dbin,dispersion,colour=gbin,group=gbin))+scale_y_log10()
ggplot(result)+geom_line(aes(dbin,exposure,colour=gbin,group=gbin))
ggplot(result)+geom_point(aes(count,dispersion))+scale_y_log10()+scale_x_log10()
ggplot(result)+geom_point(aes(mu,dispersion))+scale_y_log10()+scale_x_log10()



### group by ligation efficiency AND log10-distance bin with a single exposure
data = list( ncounts=gcounts[,.N], count=gcounts[,count], log_mean=gcounts[,log_mean],
             ntypes=gcounts[,nlevels(type)], type=gcounts[,unclass(type)],
             ngroups=gcounts[,nlevels(mixbin)], dispgroup=gcounts[,unclass(mixbin)])
op9=optimizing(stanmodel, data=data, as_vector=F, hessian=F, iter=10000, verbose=T)
op9$par$exposure
result=cbind(gcounts[,.(count=mean(count),mu=mean(exp(log_mean)),ratio=mean(count/exp(log_mean))),keyby=c("mixbin","gbin","dbin")],dispersion=op9$par$dispersion)
ggplot(result)+geom_point(aes(count,dispersion))+scale_y_log10()+scale_x_log10()
ggplot(result)+geom_point(aes(mu,dispersion))+scale_y_log10()+scale_x_log10()
ggplot(result)+geom_point(aes(ratio,dispersion))+scale_y_log10()+scale_x_log10()


### use log_dispersion = b*exp(a*log_mean) with b>0
stanmodel2 = stan_model("dispersion_diagnostics2.stan")
data = list( ncounts=gcounts[,.N], count=gcounts[,count], log_mean=gcounts[,log_mean],
             ntypes=gcounts[,nlevels(type)], type=gcounts[,unclass(type)],
             indep=gcounts[,log_mean])
op10=optimizing(stanmodel2, data=data, as_vector=F, hessian=F, iter=10000, verbose=T, init=list(exposures=array(0,dim=1),a=1,b=10))
op10$par[c("exposure","a","b")]
result=cbind(gcounts[,.(count=mean(count),mu=mean(exp(log_mean))),keyby=c("mixbin","gbin","dbin")],dispersion=op9$par$dispersion)
ggplot(result)+geom_point(aes(mu,dispersion))+scale_y_log10()+scale_x_log10()+stat_function(fun=function(x){op10$par$b*(x^op10$par$a)})



