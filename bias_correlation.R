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
data=data[2:2014]
pairs(data)


fit_nu=lm(log_nu ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
fit_delta=lm(log_delta ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
fit_npd=lm(npd ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
fit_nmd=lm(nmd ~ m40sym+m40u+m40d+gc40sym+gc40u+gc40d+gc10sym+gc10u+gc10d, data=data)
