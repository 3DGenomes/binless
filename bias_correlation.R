library(csnorm)
library(data.table)
library(ggplot2)

load("data/caulo_NcoI_all_csnorm_optimized.RData")

#produce mappability info at cut site locations
get_mappability = function(positions, bases.up=10, bases.dn=10) {
  map=fread("~/simulations/caulobacter/data/caulobacter_mappability/MQscore.dat")$V1
  sums=cumsum(map)
  ends=pmin(length(map),positions+bases.dn)
  begins=pmax(1,positions-bases.up-1)
  return((sums[ends]-sums[begins])/(bases.up+bases.dn+1))
}

map=get_mappability(cs@biases[,pos])

get_gc_content = function(positions, bases.up=10, bases.dn=10) {
  gc=readLines("~/simulations/caulobacter/data/Caulobacter_crescentus_na1000.GCA_000022005.1.25.dna.genome.fa")
  gc=do.call(c, lapply(gc[2:length(gc)],function(x){substring(x,1:nchar(x),1:nchar(x))}))
  sums=cumsum(gc %in% c("G","C"))
  ends=pmin(length(map),positions+bases.dn)
  begins=pmax(1,positions-bases.up-1)
  return((sums[ends]-sums[begins])/(bases.up+bases.dn+1))
}

gc=get_gc_content(cs@biases[,pos])

get_count_sums = function(positions, counts) {
  cs1=counts[,sum(contact.close+contact.down+contact.far+contact.up),by=pos1][,.(pos=pos1,N=V1)]
  cs2=counts[,sum(contact.close+contact.down+contact.far+contact.up),by=pos2][,.(pos=pos2,N=V1)]
  setkey(cs1,pos)
  setkey(cs2,pos)
  pos=data.table(pos=positions, key="pos")
  sums=rbind(cs1,cs2)[,sum(N),keyby=pos][pos]
  return(sums[,V1])
}

get_count_diffs = function(positions, counts) {
  cs1=counts[,sum(-contact.close-contact.down+contact.far+contact.up),by=pos1][,.(pos=pos1,N=V1)]
  cs2=counts[,sum(+contact.close-contact.down-contact.far+contact.up),by=pos2][,.(pos=pos2,N=V1)]
  setkey(cs1,pos)
  setkey(cs2,pos)
  pos=data.table(pos=positions, key="pos")
  sums=rbind(cs1,cs2)[,sum(N),keyby=pos][pos]
  return(sums[,V1])
}

get_left_right_count_sums = function(positions, counts) {
  cs1=counts[,.(R=sum(contact.close+contact.down),L=sum(contact.far+contact.up)),by=pos1][,.(pos=pos1,L,R)]
  cs2=counts[,.(R=sum(contact.far+contact.down),L=sum(contact.close+contact.up)),by=pos2][,.(pos=pos2,L,R)]
  pos=data.table(pos=positions, key="pos")
  sums=rbind(cs1,cs2)[,.(L=sum(L),R=sum(R)),keyby=pos][pos]
  return(sums)
}


get_nu_estimates = function(positions, counts) {
  S=length(positions)
  dt=data.table(pos=positions,sums=get_count_sums(positions, counts))
  eC=dt[sums>0,sum(log(sums))/S-log(4*S)]
  dt[,log_nu:=log(sums)-log(4*S)-eC]
  min_log_nu=dt[!is.infinite(log_nu),min(log_nu)]
  dt[is.infinite(log_nu),log_nu:=min_log_nu/10]
  return(list(eC=eC, log_nu=dt[,log_nu]))
}


get_nu_delta_estimates = function(positions, counts) {
  S=length(positions)
  sums=get_left_right_count_sums(positions, counts)
  sums[L==0,L:=1]
  sums[R==0,R:=1]
  sums[,c("log_delta","log_nu"):=list((log(L)-log(R))/2,(log(L)+log(R))/2)]
  sums[,c("log_delta","log_nu"):=list(log_delta-sum(log_delta)/S,log_nu-sum(log_nu)/S)]
  eC=sums[,mean(c(log(L/(2*S))-log_nu-log_delta,log(R/(2*S))-log_nu+log_delta))]
  return(list(eC=eC, log_nu=sums[,log_nu], log_delta=sums[,log_delta]))
}

get_nu_delta_estimates2 = function(positions, counts) {
  S=length(positions)
  sums=get_left_right_count_sums(positions, counts)
  sums[L==0,L:=1]
  sums[R==0,R:=1]
  sums[,c("log_delta","nu"):=list((L-R)/(L+R),(L+R)/2)]
  sums[,nu:=nu/mean(nu)]
  sums[,log_nu:=nu-1]
  eC=sums[,mean(c(log(L/(2*S))-log_nu-log_delta,log(R/(2*S))-log_nu+log_delta))]
  return(list(eC=eC, log_nu=sums[,log_nu], log_delta=sums[,log_delta]))
}


get_count_ratios = function(positions, counts) {
  cs1=counts[,sum(contact.far+contact.up)/sum(contact.close+contact.down),by=pos1][,.(pos=pos1,N=V1)]
  cs2=counts[,sum(contact.close+contact.up)/sum(contact.down+contact.far),by=pos2][,.(pos=pos2,N=V1)]
  setkey(cs1,pos)
  setkey(cs2,pos)
  pos=data.table(pos=positions, key="pos")
  sums=rbind(cs1,cs2)[,sum(N),keyby=pos][pos]
  return(sums[,V1])
}

get_bias_sums = function(biases) {
  return(biases[,dangling.L+dangling.R+rejoined])
}

get_de_diff = function(biases) {
  return(biases[,dangling.L-dangling.R])
}


data=data.table(pos=cs@biases[,pos], log_nu=cs@par$log_nu, log_delta=cs@par$log_delta)
data[,npd:=log_nu+log_delta]
data[,nmd:=log_nu-log_delta]
data[,m40sym:=1-10**(-0.1*get_mappability(pos,bases.up=20,bases.dn=20))]
data[,m40u:=1-10**(-0.1*get_mappability(pos,bases.up=40,bases.dn=0))]
data[,m40d:=1-10**(-0.1*get_mappability(pos,bases.up=0,bases.dn=40))]
data[,gc40sym:=get_gc_content(pos,bases.up=20,bases.dn=20)]
data[,gc40u:=get_gc_content(pos,bases.up=40,bases.dn=0)]
data[,gc40d:=get_gc_content(pos,bases.up=0,bases.dn=40)]
data[,gc10sym:=get_gc_content(pos,bases.up=5,bases.dn=5)]
data[,gc10u:=get_gc_content(pos,bases.up=10,bases.dn=0)]
data[,gc10d:=get_gc_content(pos,bases.up=0,bases.dn=10)]
data[,count.sums:=get_count_sums(pos,cs@counts)]
data[,count.diffs:=get_count_diffs(pos,cs@counts)]
data[,count.ratios:=get_count_ratios(pos,cs@counts)]
data[,bias.sums:=get_bias_sums(cs@biases)]
data[,de.diffs:=get_de_diff(cs@biases)]
#data=data[2:2014]
pairs(data=data, ~log_nu+log_delta+npd+nmd+count.sums+bias.sums+count.diffs+de.diffs)


fit_nu=lm(log_nu ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
fit_delta=lm(log_delta ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
fit_npd=lm(npd ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
fit_nmd=lm(nmd ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)

summary(lm(log_nu~count.sums+bias.sums+count.diffs+de.diffs,data=data))
summary(lm(log_delta~count.sums+bias.sums+count.diffs+de.diffs,data=data))
summary(lm(log_delta~count.diffs/count.sums,data=data))
summary(lm(log_delta~I(count.diffs/count.sums),data=data))
summary(lm(log_delta~count.diffs,data=data))
pairs(~I(log_delta-log_nu)+count.diffs+I(count.diffs/count.sums)+log(count.ratios),data=data[!is.infinite(count.ratios)])

data[,nu:=exp(log_nu)]
data[,delta:=exp(log_delta)]
data[,csum.std:=count.sums/mean(count.sums)]
data[,bsum.std:=bias.sums/mean(bias.sums)]
data[,cdiff.std:=count.diffs/count.sums]

#nu vs count sum
pairs(data=data, ~nu+log_nu+count.sums+csum.std+log(csum.std))
ggplot(data)+geom_point(aes(log_nu,csum.std))+xlim(-2,2)+ylim(-2,2)
summary(lm(log_nu~csum.std,data=data))
ggplot(data)+geom_point(aes(nu,csum.std))+xlim(-2,2)+ylim(-2,2)
summary(lm(nu~csum.std,data=data))
ggplot(data)+geom_point(aes(log_nu,log(csum.std)))+xlim(-2,2)+ylim(-2,2)
summary(lm(log_nu~log(count.sums),data=data))
summary(lm(log_nu~log(csum.std)+0,data=data))
summary(lm(log_nu~log(csum.std)+log(bsum.std)+0,data=data[bsum.std>0]))

#delta vs count sums and diffs
pairs(data=data, ~delta+log_delta+count.diffs+cdiff.std+count.sums)


data[,log_nu_est:=get_nu_delta_estimates(pos,cs@counts)$log_nu]
data[,log_nu_est2:=get_nu_delta_estimates2(pos,cs@counts)$log_nu]
pairs(data=data[log_nu>-1.5], ~nu+log_nu+count.sums+csum.std+log(csum.std)+log_nu_est+exp(log_nu_est)+log_nu_est2)
summary(rlm(log_nu~offset(log(csum.std)),data=data))
summary(rlm(log_nu~offset(log_nu_est),data=data))
summary(rlm(log_nu~offset(log_nu_est2),data=data))

data[,log_delta_est:=get_nu_delta_estimates(pos,cs@counts)$log_delta]
data[,log_delta_est2:=get_nu_delta_estimates2(pos,cs@counts)$log_delta]
pairs(data=data[delta<2], ~delta+log_delta+count.diffs+cdiff.std+count.sums+log_delta_est+exp(log_delta_est)+log_delta_est2)
summary(rlm(log_delta~offset(cdiff.std),data=data))
summary(rlm(log_delta~offset(log_delta_est),data=data))
summary(rlm(log_delta~offset(log_delta_est2),data=data))

sm=stan_model("guess_aggreg.stan")
bf_per_kb=0.25
sums=get_left_right_count_sums(positions, counts)
op=optimizing(sm, data=list(Krow=round(cs@biases[,(max(pos)-min(pos))/1000*bf_per_kb]), S=cs@biases[,.N],
                            cutsites=cs@biases[,pos], rejoined=cs@biases[,rejoined],
                            danglingL=cs@biases[,dangling.L], danglingR=cs@biases[,dangling.R],
                            counts_sum_left=sums[,L], counts_sum_right=sums[,R]),
              as_vector=F, hessian=F, iter=10000, verbose=T, init=0)
#
data[,log_nu_stan:=op$par$log_nu]
pairs(data=data[log_nu>-1.5&log_nu_stan<1], ~log_nu+log(csum.std)+log_nu_est2+log_nu_stan)
summary(rlm(log_nu~offset(log_nu_est2),data=data))
summary(rlm(log_nu~offset(log_nu_stan),data=data))
#
data[,log_delta_stan:=op$par$log_delta]
pairs(data=data[delta<2], ~log_delta+cdiff.std+log_delta_est2+log_delta_stan)
summary(rlm(log_nu~offset(log_delta_est2),data=data))
summary(rlm(log_nu~offset(log_delta_stan),data=data))
