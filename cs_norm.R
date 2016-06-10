library(data.table)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(shinystan)
library(mgcv)
library(scam)
library(Hmisc)
library(foreach)
library(doParallel)
library(MASS)
library(matrixStats)

setwd("/home/yannick/simulations/cs_norm")


read_tsv = function(fname, nrows=-1L, skip=0L) {
  data = fread(fname, col.names=c("id", "chr1","begin1","strand1", "length1",
                                  "re.up1","re.dn1",
                                  "chr2","begin2","strand2", "length2",
                                  "re.up2","re.dn2"), nrows=nrows, skip=skip)
  #only one chromosome for now
  stopifnot(data[,nlevels(factor(chr1))]==1)
  stopifnot(data[,nlevels(factor(chr2))]==1)
  stopifnot(data[,chr1[1]==chr2[1]])
  data[,chr1:=NULL]
  data[,chr2:=NULL]
  #
  message("put read2 always downstream of read1")
  bdata=data[begin1>begin2]
  setnames(bdata, c("begin1","strand1", "length1", "re.up1","re.dn1",
                    "begin2","strand2", "length2", "re.up2","re.dn2"),
           c("begin2","strand2", "length2", "re.up2","re.dn2",
             "begin1","strand1", "length1", "re.up1","re.dn1"))
  data = rbind(data[begin1<=begin2],bdata)
  #
  message("set rbegin and rend")
  data[strand1==0,c("rbegin1","rend1"):=list(begin1,begin1-length1+1)]
  data[strand1==1,c("rbegin1","rend1"):=list(begin1,begin1+length1-1)]
  data[strand2==0,c("rbegin2","rend2"):=list(begin2,begin2-length2+1)]
  data[strand2==1,c("rbegin2","rend2"):=list(begin2,begin2+length2-1)]
  #
  message("find each read's closest restriction site")
  data[,re.closest1:=ifelse(rbegin1-re.up1 < re.dn1-rbegin1,re.up1,re.dn1)]
  data[,re.closest2:=ifelse(rbegin2-re.up2 < re.dn2-rbegin2,re.up2,re.dn2)]
  #
  message("number it")
  rsites=data.table(pos=data[,unique(c(re.closest1,re.closest2))])
  setkey(rsites, pos)
  rsites[,idx:=.I]
  setkey(data,re.closest1)
  data=rsites[data]
  setnames(data,"pos","re.closest1")
  setnames(data,"idx","re.closest1.idx")
  setkey(data,re.closest2)
  data=rsites[data]
  setnames(data,"pos","re.closest2")
  setnames(data,"idx","re.closest2.idx")
  #
  message("sort by begins")
  setkey(data,rbegin1,rbegin2)
  return(data)
}

get_subset = function(data, b1, e1, b2=NULL, e2=NULL) {
  if (is.null(b2)) b2=b1
  if (is.null(e2)) e2=e1
  stopifnot(e1>b1)
  stopifnot(e2>b2)
  stopifnot(b2>=b1)
  return(data[begin1>=b1 & begin1+length1<=e1 & begin2>=b2 & begin2+length2<=e2])
}

categorize_by_new_type = function(sub, maxlen=600, dangling.L = c(0,4), dangling.R = c(3,-1)) {
  sub[,category:="other"]
  #careful: order is important because attribution criteria overlap
  message("Contacts")
  sub[re.closest1.idx != re.closest2.idx & re.dn1 != re.dn2,
      category:=ifelse(strand1==1 & re.dn1 - rbegin1 < maxlen & strand2==1 & re.dn2 - rbegin2 < maxlen, "contact up",
                       ifelse(strand1==1 & re.dn1 - rbegin1 < maxlen & strand2==0 & rbegin2 - re.up2 < maxlen & rbegin2-rbegin1 >= maxlen, "contact far", #necessary to not overwrite DE
                              ifelse(strand1==0 & rbegin1 - re.up1 < maxlen & strand2==1 & re.dn2 - rbegin2 < maxlen, "contact close",
                                     ifelse(strand1==0 & rbegin1 - re.up1 < maxlen & strand2==0 & rbegin2 - re.up2 < maxlen, "contact down", "other"))))]
  message("Random")
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & re.dn1 == re.dn2, category:="random"]
  message("Rejoined")
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & re.dn1 <= re.up2, category:="rejoined"]
  #sub[re.closest1.idx != re.closest2.idx & rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0, category:="random"]
  message("Self Circle")
  sub[re.dn1 == re.dn2 & strand1==0 & strand2==1 & rbegin1 - re.up1 < maxlen & re.dn2 - rbegin2 < maxlen, category:="self circle"]
  message("Dangling L/R")
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & ((rbegin1 - re.closest1) %in% dangling.L), category:="dangling L"]
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & ((rbegin2 - re.closest2) %in% dangling.R), category:="dangling R"]
  return(sub)
}

iterative_normalization = function(bdata, niterations=100, resolution=1000000, return.binned=F) {
  binned = bdata[bin1<bin2,.(bin1,bin2,N=observed)]
  binned = rbind(binned, binned[,.(bin1=bin2,bin2=bin1,N)])
  binned[,N.weighted:=N]
  #iterate
  for (i in 1:niterations) {
    binned[,b1:=sum(N.weighted),by=bin1]
    binned[,b2:=sum(N.weighted),by=bin2]
    binned[,b1:=b1/mean(b1)]
    binned[,b2:=b2/mean(b2)]
    binned[,N.weighted:=N.weighted/b1/b2]
  }
  if (return.binned == T) {
    binned[,c("b1","b2","N"):=list(NULL,NULL,NULL)]
    setnames(binned,"N.weighted","N")
    binned[,begin1:=do.call(as.integer, tstrsplit(as.character(bin1), "[[,]")[2])]
    binned[,end1:=do.call(as.integer, tstrsplit(as.character(bin1), "[],)]")[2])]
    binned[,begin2:=do.call(as.integer, tstrsplit(as.character(bin2), "[[,]")[2])]
    binned[,end2:=do.call(as.integer, tstrsplit(as.character(bin2), "[],)]")[2])]
    return(binned[begin1<begin2])
  } else {
    #get weights and report them in data
    binned[,weight:=N.weighted/N]
    setkey(binned, bin1, bin2)
    setkey(bdata, bin1, bin2)
    data[category %in% c("contact close","contact up","contact down","contact far"),weight:=binned[bdata, weight]]
    return(data)
  }
}

prepare_for_sparse_cs_norm = function(data, both=T, circularize=-1) {
  #get counts that are not other, random or self circle
  message("Sum up counts and biases")
  enrich = data[!(category %in% c("other","random", "self circle")), .N,
                keyby=c("re.closest1","re.closest2", "category")]
  enrich[,category:=gsub(" ", ".", category)]
  #
  #dangling ends: attribute to cut site
  message("Reattribute biases to proper rsites")
  enrich[category=="dangling.L", re.closest2:=re.closest1]
  enrich[category=="dangling.R", re.closest1:=re.closest2]
  #rejoined ends: attribute to rsite in the middle
  enrich[category=="rejoined", re.pos:=as.integer((re.closest1+re.closest2)/2.)]
  rsites = data.table(re.pos=enrich[,unique(c(re.closest1,re.closest2))])
  setkey(enrich,re.pos)
  setkey(rsites,re.pos)
  enrich=rsites[,.(re.pos,re.closest=re.pos)][enrich,,roll="nearest"]
  enrich[category=="rejoined",c("re.closest1","re.closest2"):=list(re.closest,re.closest)]
  enrich[,c("re.pos","re.closest"):=list(NULL,NULL)]
  #
  #get rsites that a) have counts on them and b) have at least one dangling/rejoined end.
  message("RSites list")
  rsites.count=enrich[category %in% c("contact.close", "contact.down", "contact.far", "contact.up"),
                      unique(c(re.closest1,re.closest2))]
  rsites.bias=enrich[category %in% c("dangling.L", "dangling.R", "rejoined"),
                     unique(c(re.closest1,re.closest2))]
  rsites=data.table(re.pos=intersect(rsites.count, rsites.bias))
  setkey(rsites, re.pos)
  rsites[,id:=.I]
  #
  #attribute restriction site IDs
  message("Add Rsite info")
  enrich = merge(enrich, rsites, by.x="re.closest1", by.y="re.pos")
  enrich = merge(enrich, rsites, by.x="re.closest2", by.y="re.pos", suffixes=c("1", "2"))
  #
  message("Produce bias list")
  X=dcast.data.table(enrich[id1 == id2 & category %in% c("dangling.L", "dangling.R", "rejoined"),
                            .(id=id1, pos=re.closest1, category, N)],
                     ...~category, value.var="N", fun.aggregate = sum)
  stopifnot(all(X[,id==.I])) #IDs must start at one and have no gaps, for current stan implementation
  setkey(X,id)
  #
  message("Produce count list")
  Y=dcast.data.table(enrich[id1 != id2 
                            & category %in% c("contact.close", "contact.down", "contact.far", "contact.up"),
                            .(id1,id2,pos1=re.closest1,pos2=re.closest2,category,N)],
                     ...~category, value.var="N", fill=0)
  setkey(Y,id1,id2)
  Y[,distance:=pos2-pos1]
  if (circularize>0) Y[,distance:=pmin(distance,circularize-distance+1)]
  #
  if (both == T) {
    message("Merge them")
    Ymelt = melt(Y, id.vars = c("id1","pos1", "id2", "pos2", "distance"), variable.name = "category", value.name = "N")
    XY = merge(Ymelt, X, by.x=c("id1","pos1"), by.y=c("id","pos"))
    XY = merge(XY, X, by.x=c("id2","pos2"), by.y=c("id","pos"), suffixes=c(".1", ".2"))
    XY = merge(XY, rsites, by.x=c("id1","pos1"), by.y=c("id","re.pos"))
    XY = merge(XY, rsites, by.x=c("id2","pos2"), by.y=c("id","re.pos"), suffixes = c("1","2"))
    setkey(XY, id1, id2)
  } else {
    XY=NULL
  }
  return(list(biases=X, counts=Y, both=XY))
}

bin_for_mean_field = function(biases, counts, distance_bins_per_decade=100, circularize=-1) {
  stopifnot(counts[id1>=id2,.N]==0)
  mcounts=melt(counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0]
  ### accumulate counts
  ci=mcounts[,.(id=id1,count,category)][,.N,by=c("id","count","category")]
  cj=mcounts[,.(id=id2,count,category)][,.N,by=c("id","count","category")]
  nsites=biases[,.N]
  ### make histograms for biases
  Nkl=dcast(rbind(ci[category=="contact.up",.(id,count,category="Ni.up",N)],
                  ci[category=="contact.far",.(id,count,category="Ni.far",N)],
                  cj[category=="contact.up",.(id,count,category="Nj.up",N)],
                  cj[category=="contact.close",.(id,count,category="Nj.close",N)]),
            ...~category, value.var="N", fill=0)[,.(id,count,N=Ni.far+Ni.up+Nj.up+Nj.close)]
  Nkl=rbind(Nkl,Nkl[,.(count=0,N=2*nsites-sum(N)),by=id]) #each rsite is counted twice
  setkey(Nkl,N,id,count)
  Nkl=Nkl[N>0]
  Nkr=dcast(rbind(ci[category=="contact.close",.(id,count,category="Ni.close",N)],
                  ci[category=="contact.down",.(id,count,category="Ni.down",N)],
                  cj[category=="contact.far",.(id,count,category="Nj.far",N)],
                  cj[category=="contact.down",.(id,count,category="Nj.down",N)]),
            ...~category, value.var="N", fill=0)[,.(id,count,N=Ni.close+Ni.down+Nj.far+Nj.down)]
  Nkr=rbind(Nkr,Nkr[,.(count=0,N=2*nsites-sum(N)),by=id]) #each rsite is counted twice
  setkey(Nkr,N,id,count)
  Nkr=Nkr[N>0]
  ### make histogram for distance
  #make distance bins and their factor
  stepsz=1/distance_bins_per_decade
  if (circularize>0) {
    dmax=biases[,max(pmin(pos,circularize-pos+1))]
    dmin=biases[,min(pmin(pos,circularize-pos+1))]
    dbins=10**seq(0,log10(dmax-dmin)+stepsz,stepsz)
  } else {
    dbins=10**seq(0,biases[,log10(max(pos)-min(pos))]+stepsz,stepsz)
  }
  mcounts[,distance:=abs(pos1-pos2)]
  mcounts[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #Count positive counts in these bins
  Nkd = mcounts[,.N,keyby=c("bdist","count")]
  Nkd[,mdist:=sqrt(dbins[unclass(bdist)+1]*dbins[unclass(bdist)])]
  #Count the number of crossings per distance bin
  positions=biases[,pos]
  npos=length(positions)
  ncrossings <- rowSums(sapply(1:(npos-1), #this loop is still 5x faster in python
                               function(i){hist(positions[(i+1):npos]-positions[i],breaks=dbins,plot=F,
                                                right=F, include.lowest=T)$counts}
  ))
  #deduce zero counts
  Nkz = data.table(bdist=Nkd[,ordered(levels(bdist), levels(bdist))], ncrossings=ncrossings, key="bdist")
  Nkz = Nkd[,.(nnz=sum(N)),by=bdist][Nkz[ncrossings>0]]
  Nkz[,mdist:=sqrt(dbins[unclass(bdist)+1]*dbins[unclass(bdist)])]
  Nkz[is.na(nnz),nnz:=0]
  Nkd = rbind(Nkd, Nkz[,.(bdist,mdist,count=0,N=4*ncrossings-nnz)]) # one crossing for each of 4 count types
  setkey(Nkd,N,bdist,count)
  Nkd=Nkd[N>0]
  #plot decay
  #ggplot(Nkd[,.(decay=sum(N*count)/sum(N)),by=bdist])+geom_point(aes(bdist,decay))+scale_y_log10()+
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  stopifnot(Nkr[,all(N>=0)],Nkl[,all(N>=0)],Nkd[,all(N>=0)])
  return(list(Nkd=Nkd, Nkr=Nkr, Nkl=Nkl))
}

dset_statistics = function(biases,counts){
  message("Mean counts")
  message("   Rejoined  : ", biases[,mean(rejoined)])
  message("   Dangling L: ", biases[,mean(dangling.L)])
  message("   Dangling R: ", biases[,mean(dangling.R)])
  message("   C. close  : ", counts[,mean(contact.close)])
  message("   C. far    : ", counts[,mean(contact.far)])
  message("   C. up     : ", counts[,mean(contact.up)])
  message("   C. down   : ", counts[,mean(contact.down)])
  message("Median counts")
  message("   Rejoined  : ", biases[,median(rejoined)])
  message("   Dangling L: ", biases[,median(dangling.L)])
  message("   Dangling R: ", biases[,median(dangling.R)])
  message("   C. close  : ", counts[,median(contact.close)])
  message("   C. far    : ", counts[,median(contact.far)])
  message("   C. up     : ", counts[,median(contact.up)])
  message("   C. down   : ", counts[,median(contact.down)])
  message("Percent of zero counts")
  message("   Rejoined  : ", biases[rejoined==0,.N]/biases[,.N]*100)
  message("   Dangling L: ", biases[dangling.L==0,.N]/biases[,.N]*100)
  message("   Dangling R: ", biases[dangling.R==0,.N]/biases[,.N]*100)
  message("   C. close  : ", counts[contact.close==0,.N]/counts[,.N]*100)
  message("   C. far    : ", counts[contact.far==0,.N]/counts[,.N]*100)
  message("   C. up     : ", counts[contact.up==0,.N]/counts[,.N]*100)
  message("   C. down   : ", counts[contact.down==0,.N]/counts[,.N]*100)
}

generate_fake_dataset = function(num_rsites=10, genome_size=10000, eC=.1, eRJ=.4, eDE=.8, alpha=10) {
  #place rsites
  biases=data.table(id=seq(num_rsites),
                    pos=cumsum(rmultinom(n=1, size=genome_size, prob=rep(1,num_rsites+1)))[1:num_rsites])
  setkey(biases,id)
  #build biases
  biases[,true_log_nu:=sin(pos/1000)+(pos-genome_size/2)/genome_size]
  biases[,true_log_nu:=true_log_nu-mean(true_log_nu)]
  #ggplot(biases,aes(pos,true_log_nu))+geom_point()+geom_line()
  biases[,true_log_delta:=-(pos-genome_size/2)/genome_size+sin(pos/2100)]
  biases[,true_log_delta:=true_log_delta-mean(true_log_delta)]
  #ggplot(biases,aes(pos,true_log_delta))+geom_point()+geom_line()
  biases[,true_log_mean_RJ:=eRJ+true_log_nu]
  biases[,true_log_mean_DL:=eDE+true_log_nu+true_log_delta]
  biases[,true_log_mean_DR:=eDE+true_log_nu-true_log_delta]
  #draw dangling/rejoined
  biases[,dangling.L:=rnbinom(.N, mu=exp(true_log_mean_DL), size=alpha)]
  biases[,dangling.R:=rnbinom(.N, mu=exp(true_log_mean_DR), size=alpha)]
  biases[,rejoined:=rnbinom(.N, mu=exp(true_log_mean_RJ), size=alpha)]
  #report rsites in counts
  counts=CJ(biases[,paste(id,pos,true_log_nu,true_log_delta)],biases[,paste(id,pos,true_log_nu,true_log_delta)])
  counts[,c("id1","pos1","true_log_nu1","true_log_delta1"):=tstrsplit(V1, " ")]
  counts[,c("id2","pos2","true_log_nu2","true_log_delta2"):=tstrsplit(V2, " ")]
  counts[,c("id1","id2","pos1","pos2","V1","V2"):=list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
  counts[,c("true_log_nu1","true_log_nu2","true_log_delta1","true_log_delta2"):=list(
    as.numeric(true_log_nu1), as.numeric(true_log_nu2), as.numeric(true_log_delta1), as.numeric(true_log_delta2))]
  counts=counts[pos1<pos2]
  setkey(counts, id1, id2)
  #build decay
  counts[,distance:=pos2-pos1]
  counts[,true_log_decay:=5*dnorm(log(distance), mean=log(min(distance)), sd=log(max(distance))/10)-0.2*log(distance)]
  counts[,true_log_decay:=true_log_decay-mean(true_log_decay)]
  #ggplot(counts)+geom_point(aes(distance,exp(true_log_decay)))+scale_x_log10()+scale_y_log10()
  counts[,base_count:=eC+true_log_decay+true_log_nu1+true_log_nu2]
  counts[,true_log_mean_cclose:=base_count-true_log_delta1+true_log_delta2]
  counts[,true_log_mean_cfar:=base_count+true_log_delta1-true_log_delta2]
  counts[,true_log_mean_cup:=base_count+true_log_delta1+true_log_delta2]
  counts[,true_log_mean_cdown:=base_count-true_log_delta1-true_log_delta2]
  #draw counts
  counts[,contact.close:=rnbinom(.N, mu=exp(true_log_mean_cclose), size=alpha)]
  counts[,contact.far:=rnbinom(.N, mu=exp(true_log_mean_cfar), size=alpha)]
  counts[,contact.up:=rnbinom(.N, mu=exp(true_log_mean_cup), size=alpha)]
  counts[,contact.down:=rnbinom(.N, mu=exp(true_log_mean_cdown), size=alpha)]
  #ggplot(dset$counts)+geom_point(aes(distance,contact.close/exp(true_log_mean_cclose-true_log_decay)), alpha=0.1)+
  #  geom_line(aes(distance,exp(true_log_decay)))+scale_x_log10()+scale_y_log10()
  #statistics
  dset_statistics(biases,counts)
  return(list(biases=biases, counts=counts))
}

fill_zeros = function(counts,biases,biases2=NULL) {
  if (is.null(biases2)) biases2=biases
  newcounts=CJ(biases[,paste(id,pos)],biases2[,paste(id,pos)])
  newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
  newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
  newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
  newcounts=newcounts[pos1<pos2]
  setkey(newcounts, id1, id2, pos1, pos2)
  setkey(counts, id1, id2, pos1, pos2)
  newcounts=counts[newcounts]
  newcounts[is.na(contact.close),contact.close:=0]
  newcounts[is.na(contact.far),contact.far:=0]
  newcounts[is.na(contact.up),contact.up:=0]
  newcounts[is.na(contact.down),contact.down:=0]
  return(newcounts)
}

stan_matrix_to_datatable = function(opt, x) {
  vals=data.table(opt)
  vals[,x:=x]
  melt(data.table(vals), id.vars="x")
}

optimize_all_meanfield = function(model, biases, counts, meanfield, maxcount, bf_per_kb=1, bf_per_decade=5,
                                  iter=10000, verbose=T, init=0, ...) {
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  cclose=counts[contact.close>maxcount,.(id1,id2,distance,count=contact.close)]
  cfar=counts[contact.far>maxcount,.(id1,id2,distance,count=contact.far)]
  cup=counts[contact.up>maxcount,.(id1,id2,distance,count=contact.up)]
  cdown=counts[contact.down>maxcount,.(id1,id2,distance,count=contact.down)]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               dmin=dmin, dmax=dmax, Kdiag=Kdiag,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               Nl=mf$Nkl[,.N], Nkl_count=mf$Nkl[,count], Nkl_cidx=mf$Nkl[,id], Nkl_N=mf$Nkl[,N], Nkl_levels=mf$Nkl[,sum(diff(N)!=0)+1],
               Nr=mf$Nkr[,.N], Nkr_count=mf$Nkr[,count], Nkr_cidx=mf$Nkr[,id], Nkr_N=mf$Nkr[,N], Nkr_levels=mf$Nkr[,sum(diff(N)!=0)+1],
               Nd=mf$Nkd[,.N], Nkd_count=mf$Nkd[,count], Nkd_d=mf$Nkd[,mdist], Nkd_N=mf$Nkd[,N], Nkd_levels=mf$Nkd[,sum(diff(N)!=0)+1])
  if (data$Nl==0) data$Nkl_levels=0
  if (data$Nr==0) data$Nkr_levels=0
  if (data$Nd==0) data$Nkd_levels=0
  message("Mean field optimization")
  message("Krow        : ", Krow)
  message("Kdiag       : ", Kdiag)
  message("Biases      : ", biases[,.N])
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  message("Left counts : ", mf$Nkl[,.N], " (", data$Nkl_levels, " levels)")
  message("Right counts: ", mf$Nkr[,.N], " (", data$Nkr_levels, " levels)")
  message("Decay counts: ", mf$Nkd[,.N], " (", data$Nkd_levels, " levels)")
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

predict_all_meanfield = function(model, biases, counts, meanfield, opt, bf_per_decade=5, verbose=T) {
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag_centered=opt$par$beta_diag_centered)
  message("Mean field prediction")
  message("Kdiag       : ", Kdiag)
  message("Close counts: ", cclose[,.N])
  message("Far counts  : ", cfar[,.N])
  message("Up counts   : ", cup[,.N])
  message("Down counts : ", cdown[,.N])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
}

get_binned_matrices = function(model, biases, counts, meanfield, opt, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL, bf_per_decade=5, verbose=F) {
  stopifnot(counts[distance!=pos2-pos1,.N]==0) #need to implement circular genomes
  csub=copy(counts) #need to implement taking only needed part of matrix, and reporting log_nu and log_delta appropriately
  bsub=copy(biases)
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  Kdiag=round(csub[,(log10(dmax)-log10(dmin))*bf_per_decade])
  npoints=100*Kdiag #evaluate spline with 100 equidistant points per basis function
  #bin existing counts and biases
  if (is.null(b1)) b1=bsub[,min(pos)]-1
  if (is.null(b2)) b2=bsub[,min(pos)]-1
  if (is.null(e1)) e1=bsub[,max(pos)]+1
  if (is.null(e2)) e2=bsub[,max(pos)]+1
  bins1=seq(b1,e1+resolution,resolution)
  bins2=seq(b2,e2+resolution,resolution)
  csub[,c("bin1","bin2"):=list(cut(pos1, bins1, ordered_result=T, right=F, include.lowest=T,dig.lab=5),
                                 cut(pos2, bins2, ordered_result=T, right=F, include.lowest=T,dig.lab=5))]
  csub[,log_mean_cup:=opt$pred$log_mean_cup]
  csub[,log_mean_cdown:=opt$pred$log_mean_cdown]
  csub[,log_mean_cclose:=opt$pred$log_mean_cclose]
  csub[,log_mean_cfar:=opt$pred$log_mean_cfar]
  bsub[,log_nu:=opt$par$log_nu]
  bsub[,log_delta:=opt$par$log_delta]
  bsub[,c("bin1","bin2"):=list(cut(pos, bins1, ordered_result=T, right=F, include.lowest=T,dig.lab=5),
                                 cut(pos, bins2, ordered_result=T, right=F, include.lowest=T,dig.lab=5))]
  bsub[,c("ibin1","ibin2"):=list(as.integer(bin1),as.integer(bin2))]
  #run computation of mean
  cclose=csub[,.(id1,id2,bin1,bin2,distance,count=contact.close,logmean=log_mean_cclose)]
  cfar=csub[,.(id1,id2,bin1,bin2,distance,count=contact.far,logmean=log_mean_cfar)]
  cup=csub[,.(id1,id2,bin1,bin2,distance,count=contact.up,logmean=log_mean_cup)]
  cdown=csub[,.(id1,id2,bin1,bin2,distance,count=contact.down,logmean=log_mean_cdown)]
  csub=rbind(cclose,cfar,cup,cdown)[!(is.na(bin1)|is.na(bin2)|count==0)]
  data = list( Kdiag=Kdiag, 
               S1=bsub[!is.na(bin1),.N], S2=bsub[!is.na(bin2),.N], 
               cutsites1=bsub[!is.na(bin1),pos], cutsites2=bsub[!is.na(bin2),pos],
               dmin=dmin, dmax=dmax, npoints=npoints,
               N=csub[,.N],   counts=csub[,count], cdist=csub[,distance], cmean=csub[,exp(logmean)],
               eC=opt$par$eC,
               log_nu1=bsub[!is.na(bin1),log_nu], log_nu2=bsub[!is.na(bin2),log_nu],
               log_delta1=bsub[!is.na(bin1),log_delta], log_delta2=bsub[!is.na(bin2),log_delta],
               beta_diag_centered=opt$par$beta_diag_centered,
               B1=csub[,nlevels(bin1)], B2=csub[,nlevels(bin2)],
               cbins1=csub[,as.integer(bin1)], cbins2=csub[,as.integer(bin2)],
               bbins1=bsub[!is.na(bin1),ibin1], bbins2=bsub[!is.na(bin2),ibin2])
  binned=optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
  #format output
  so = data.table(melt(binned$par$observed))
  setnames(so,"value","observed")
  so[,expected:=melt(binned$par$expected)$value]
  so[,ncounts:=melt(binned$par$ncounts)$value]
  so[,normalized:=melt(binned$par$normalized)$value]
  so=so[ncounts>0]
  setkey(so,Var1)
  so=bsub[,.(bin1=bin1[1]),keyby=ibin1][so]
  setkey(so,Var2)
  so=bsub[,.(bin2=bin2[1]),keyby=ibin2][so]
  so[,c("ibin1","ibin2"):=list(NULL,NULL)]
  so[,lFC:=log2(observed/expected)]
  return(list(distance=exp(binned$par$log_dist), decay=exp(binned$par$log_decay), mat=so))
}

get_dispersions = function(model, binned, iter=10000, verbose=T) {
  data=list(B=binned[,.N],observed=binned[,observed],expected=binned[,expected],ncounts=binned[,ncounts])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
}

#compute p(gamma2>gamma1) = \int_{0}^{+infty} dx p_gamma2(x) \int_{0}^{x} dy p_gamma1(y)
compute_gamma_overlap = function(alpha1,beta1,alpha2,beta2, bounds=5, ncores=1) {
  registerDoParallel(cores=ncores)
  foreach (a1=alpha1, b1=beta1, a2=alpha2, b2=beta2, .packages="stats", .combine=c) %dopar% {
    #infinite integral does not work very well so truncate outer integral
    mu1=a1/b1
    mu2=a2/b2
    sd1=sqrt(a1)/b1
    sd2=sqrt(a2)/b2
    xmin = max(0,min(mu1-bounds*sd1,mu2-bounds*sd2))
    xmax = max(mu1+bounds*sd1,mu2+bounds*sd2)
    a=integrate(function(x){exp(dgamma(x,a2,rate=b2,log=T)+pgamma(x,a1,rate=b1,log.p=T))},xmin,xmax)
    if (a$abs.error<=0) {NA} else {a$value}
  }
}

detect_interactions = function(binned, dispersions, threshold=0.95, ncores=1){
  #report gamma parameters
  mat=copy(binned)
  mat[,c("alpha1","beta1"):=list(dispersions,dispersions/expected)]
  mat[,c("alpha2","beta2"):=list(alpha1+observed,beta1+1)]
  mat[,prob.observed.gt.expected:=compute_gamma_overlap(alpha1,beta1,alpha2,beta2,ncores=ncores)]
  mat[,is.interaction:=prob.observed.gt.expected>threshold | 1-prob.observed.gt.expected>threshold]
  #write begins/ends
  bin1.begin=mat[,bin1]
  bin1.end=mat[,bin1]
  bin2.begin=mat[,bin2]
  bin2.end=mat[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[[,]")[2][[1]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[[,)]")[2][[1]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[[,]")[2][[1]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[[,)]")[2][[1]]
  mat[,begin1:=as.integer(as.character(bin1.begin))]
  mat[,end1:=as.integer(as.character(bin1.end))]
  mat[,begin2:=as.integer(as.character(bin2.begin))]
  mat[,end2:=as.integer(as.character(bin2.end))]
  return(mat)
}

#estimates the values of the count or dispersion required to cross a given threshold
thresholds_estimator = function(observed, expected, dispersion, threshold=0.95, counts.range=c(expected/100,expected*100), disp.range=c(0,100*dispersion), compute=T) {
  #plot gamma distributions
  a1=dispersion
  b1=dispersion/expected
  a2=a1+observed
  b2=b1+1
  mu1=a1/b1
  mu2=a2/b2
  sd1=sqrt(a1)/b1
  sd2=sqrt(a2)/b2
  message("prior:     mean ",mu1, " sd ",sd1)
  message("posterior: mean ",mu2, " sd ",sd2)
  xmin = min(c(0,counts.range))
  xmax = max(counts.range)
  message("evaluation bounds: [",xmin,",",xmax,"]")
  p1=ggplot(data.table(x=c(xmin,xmax)),aes(x))+
    stat_function(fun=function(x){dgamma(x,a1,rate=b1)},aes(colour="prior"))+
    stat_function(fun=function(x){dgamma(x,a2,rate=b2)},aes(colour="posterior"))+
    scale_colour_manual("Gamma", values = c("blue", "red"))
  #helper function  
  gammaQuery=function(mu,theta,count){
    alpha1=theta
    alpha2=theta+count
    beta1=theta/mu
    beta2=beta1+1
    return(function(x){exp(dgamma(x,alpha2,rate=beta2,log=T)+pgamma(x,alpha1,rate=beta1,log.p=T))})
  }
  #query counts
  igammaQuery_c=function(x){integrate(gammaQuery(expected,dispersion,count=x),xmin,xmax)$value-threshold}
  p2=ggplot(data.table(x=c(xmin,xmax),y=c(0,1)),aes(x))+
    stat_function(fun=Vectorize(function(x){igammaQuery_c(x)+threshold}))+
    geom_hline(aes(yintercept=threshold))+ geom_hline(aes(yintercept=1-threshold))+
    labs(x="observed count", y="P(observed>expected)")
  if (compute){
    clo = uniroot(function(x){igammaQuery_c(x)+2*threshold-1},c(xmin,expected))$root
    chi = uniroot(igammaQuery_c,c(expected,xmax))$root
    message("At a threshold of ",threshold," significant counts are lower than ",clo," and higher than ",chi)
    p2=p2+geom_vline(aes(xintercept=clo))+ geom_vline(aes(xintercept=chi))
  }
  #query dispersion
  igammaQuery_theta=function(x){integrate(gammaQuery(expected,x,observed),xmin,xmax)$value-threshold}
  p3=ggplot(data.table(x=disp.range,y=c(0,1)),aes(x))+
    stat_function(fun=Vectorize(function(x){igammaQuery_theta(x)+threshold}))+
    geom_hline(aes(yintercept=threshold))+geom_hline(aes(yintercept=1-threshold))+
    labs(x="dispersion", y="P(observed>expected)")+scale_x_log10()
  if (compute) {
    dstar = uniroot(igammaQuery_theta,disp.range)$root
    p3=p3+geom_vline(aes(xintercept=dstar))
    if (observed > expected) {
      message("At a threshold of ",threshold," the dispersion must be larger than ",dstar," to reach significance")
    } else {
      message("At a threshold of ",threshold," the dispersion must be smaller than ",dstar," to reach significance")
    }
  }
  return(list(p1,p2,p3))
}

read_and_prepare = function(infile, outprefix, skip=0L, both=T, distance_bins_per_decade=100, circularize=-1) {
  message("*** READ")
  data=read_tsv(infile, skip=skip)
  message("*** CATEGORIZE")
  data = categorize_by_new_type(data)
  message("*** BIASES AND COUNTS")
  cs_data = prepare_for_sparse_cs_norm(data, both=both, circularize=circularize)
  dset_statistics(cs_data$biases,cs_data$counts)
  message("*** MEANFIELD")
  cs_data$meanfield=bin_for_mean_field(cs_data$biases, cs_data$counts, distance_bins_per_decade=distance_bins_per_decade)
  message("*** WRITE")
  write.table(cs_data$biases, file = paste0(outprefix,"_biases.dat"), quote=F, row.names = F)
  write.table(cs_data$counts, file = paste0(outprefix,"_counts.dat"), quote=F, row.names = F)
  if (both==T) {
    both=cs_data$both
    save(both, file = paste0(outprefix,"_both.RData"))
  }
  meanfield=cs_data$meanfield
  save(meanfield, file = paste0(outprefix, "_meanfield_", distance_bins_per_decade, ".RData"))
  return(cs_data)
}

get_cs_subset = function(counts, biases, begin1, end1, begin2=NULL, end2=NULL, fill.zeros=T) {
  stopifnot(counts[distance!=pos2-pos1,.N]==0)
  stopifnot(end1>=begin1)
  if (is.null(begin2)) begin2=begin1
  if (is.null(end2)) end2=end1
  stopifnot(begin2>=begin1)
  #
  biases.1=biases[pos>=begin1&pos<end1]
  beginrange1=biases.1[1,id]
  endrange1=biases.1[.N,id]
  biases.1[,id:=id-beginrange1+1]
  #
  biases.2=biases[pos>=begin2&pos<end2]
  beginrange2=biases.2[1,id]
  endrange2=biases.2[.N,id]
  biases.2[,id:=id-beginrange2+1]
  #
  counts.local=counts[pos1>=begin1&pos1<end1&pos2>=begin2&pos2<end2]
  counts.local[,c("id1","id2"):=list(id1-beginrange1+1,id2-beginrange2+1)]
  #fill zeros
  if (fill.zeros==T) {
    counts.local=fill_zeros(counts=counts.local,biases=biases.1,biases2=biases.2)
    counts.local[,distance:=abs(pos2-pos1)] #not filled
  }
  return(list(counts=counts.local,biases1=biases.1,biases2=biases.2,
              beginrange1=beginrange1, endrange1=endrange1, beginrange2=beginrange2, endrange2=endrange2))
}

optimize_genomic_params = function(model, dt, lambda_nu, lambda_delta, bf_per_kb=1, iter=10000, verbose=T, init=0, ...) {
  Krow=round(dt[,(max(pos)-min(pos))/1000*bf_per_kb])
  data = list( Krow=Krow, S=dt[,.N], cutsites=dt[,pos],
               log_mean_RJ=dt[,log_mean_RJ], log_mean_DL=dt[,log_mean_DL], log_mean_DR=dt[,log_mean_DR],
               lambda_nu=lambda_nu, lambda_delta=lambda_delta)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

optimize_extradiag = function(model, biases1, biases2, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Kdiag=Kdiag, S1=biases1[,.N], S2=biases2[,.N], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               log_nu1=biases1[,log_nu], log_nu2=biases2[,log_nu], log_delta1=biases1[,log_delta], log_delta2=biases2[,log_delta])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

optimize_decay_params = function(model, dt, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  dist=dt[,(begin+end)/2]
  fij=dt[,(log_decay_close+log_decay_far+log_decay_up+log_decay_down)/4]
  ncounts=dt[,ncounts]
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=length(dist), dist=dist, fij=fij, ncounts=ncounts)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

run_split_parallel_squares = function(biases, square.size, coverage) {
  minpos=biases[,min(pos)]
  maxpos=biases[,max(pos)+1]
  true.square.size=(maxpos-minpos)/round((maxpos-minpos)/square.size)
  begins=seq(minpos,maxpos-true.square.size,by=true.square.size/coverage)
  squares=CJ(begins,begins)[V2>=V1]
  setnames(squares, c("begin1","begin2"))
  squares[,c("end1","end2"):=list(begin1+true.square.size,begin2+true.square.size)]
  diagsquares=squares[begin1==begin2]
  #ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2),alpha=0.1)+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=factor((begin2-begin1)/true.square.size)),alpha=0.1)+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin2-end1,xmax=end2-begin1,ymin=sqid,ymax=sqid+1),alpha=0.1)+scale_x_log10()
  message("Run in parallel")
  message("   square size: ",true.square.size)
  message("   coverage: ",coverage)
  message("   ",squares[,.N], " squares")
  message("   ",diagsquares[,.N], " on the diagonal")
  return(list(squares=squares,diagsquares=diagsquares,true.square.size=true.square.size))
}

run_split_parallel_biases = function(smfit, counts, biases, meanfield, begin, end, bf_per_kb, bf_per_decade, verbose, iter) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases, begin1=begin, end1=end, fill.zeros=T)
  #run fit
  op <- optimize_all_meanfield(smfit, extracted$biases1, extracted$counts, meanfield, maxcount=-1, bf_per_kb = bf_per_kb,
                               bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)
  #report data
  center=(end-begin)/2
  ret=extracted$biases1[,.(id=id+extracted$beginrange1-1,pos,weight=(center-abs(pos-center))**2)]
  ret[,c("log_mean_DL", "log_mean_DR", "log_mean_RJ"):=list(op$par$log_mean_DL, op$par$log_mean_DR, op$par$log_mean_RJ)]
  ret[,square.begin:=begin]
  ret[,c("lambda_delta","lambda_nu","nsites","alpha"):=list(op$par$lambda_delta,op$par$lambda_nu,extracted$biases1[,.N],op$par$alpha)]
}

optimize_eC = function(model, biases, counts, retlist, dmin, dmax, bf_per_decade=5, verbose=T, iter=1000) {
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round(counts[,(log10(dmax)-log10(dmin))*bf_per_decade])
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               log_nu=retlist$log_nu, log_delta=retlist$log_delta, beta_diag_centered=retlist$beta_diag_centered)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
}

run_split_parallel = function(counts, biases, square.size=100000, coverage=4, bf_per_kb=1,
                              bf_per_decade=5, distance_bins_per_decade=100, verbose = F, iter=100000, ncpus=30, homogenize=F) {
  ### build squares
  message("*** build squares")
  retval.squares = run_split_parallel_squares(biases, square.size, coverage)
  squares=retval.squares$squares
  diagsquares=retval.squares$diagsquares
  true.square.size=retval.squares$true.square.size
  ### fit genomic biases using squares close to diagonal
  message("*** fit genomic biases")
  registerDoParallel(cores=ncpus)
  smfit = stan_model(file = "cs_norm_fit.stan")
  ops.bias = foreach (begin=diagsquares[,begin1], end=diagsquares[,end1], .packages=c("data.table","rstan"), .combine=rbind) %dopar% 
    run_split_parallel_biases(smfit, counts, biases, meanfield, begin, end, bf_per_kb = bf_per_kb,
                              bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)
  ### reconstruct bias estimates: eRJ eDE log_nu and log_delta
  message("*** reconstruct genomic biases")
  info=ops.bias[,.SD[1],by=square.begin,.SDcols=c("lambda_nu","lambda_delta","nsites","alpha")]
  #ggplot(info)+geom_line(aes(square.begin,lambda_nu))+scale_y_log10()+geom_hline(yintercept=info[,exp(median(log(lambda_nu)))])
  means=ops.bias[,.(log_mean_RJ=weighted.mean(log_mean_RJ, weight),
             log_mean_DL=weighted.mean(log_mean_DL, weight),
             log_mean_DR=weighted.mean(log_mean_DR, weight),
             ncounts=.N), keyby=c("id","pos")]
  if (homogenize==T) {
    message("*** homogenize genomic biases")
    smfitgen=stan_model("cs_norm_fit_genomic_params.stan")
    op=optimize_genomic_params(smfitgen, means, lambda_nu=info[,exp(mean(log(lambda_nu)))], lambda_delta=info[,exp(mean(log(lambda_delta)))],
                               bf_per_kb=bf_per_kb, iter=iter, verbose=verbose)
    #ggplot(means)+geom_line(aes(pos,log_mean_RJ,colour="ori"))+
    #  geom_line(data=data.table(pos=means[,pos],rj=op$par$log_nu+op$par$eRJ),aes(pos,rj,colour="avg"))
    #ggplot(means)+geom_line(aes(pos,log_mean_DL,colour="ori"))+
    #  geom_line(data=data.table(pos=means[,pos],rj=op$par$log_nu+op$par$eDE+op$par$log_delta),aes(pos,rj,colour="avg"))
    retlist=op$par[c("eRJ","eDE","log_nu","log_delta")]
    setkey(biases, id, pos)
    biases.aug=copy(biases)
    biases.aug[,log_nu:=retlist$log_nu]
    biases.aug[,log_delta:=retlist$log_delta]
  } else {
    setkey(biases, id, pos)
    retlist=list(eRJ=ops.bias[,mean(log_mean_RJ)], eDE=ops.bias[,mean(log_mean_DL+log_mean_DR)/2])
    biases.aug=means[biases]
    stopifnot(!any(is.na(biases.aug)))
    biases.aug[,log_nu:=log_mean_RJ-retlist$eRJ]
    biases.aug[,log_delta:=(log_mean_DL-log_mean_DR)/2]
    retlist$log_nu=biases.aug[,log_nu]
    retlist$log_delta=biases.aug[,log_delta]
  }
  #ggplot()+geom_line(data=data.table(x=biases.aug[,pos],y=opall$par$log_delta),aes(x,y,colour="ref"))+
  #  geom_line(data=data.table(x=biases.aug[,pos],y=bias.params$log_delta), aes(x,y,colour="split"))
  #
  ### fit remaining data
  message("*** fit diagonal decay")
  smfitex=stan_model("cs_norm_fit_extradiag.stan")
  dmax=max(counts[,max(distance)],meanfield$Nkd[,max(mdist)])+0.01
  dmin=min(counts[,min(distance)],meanfield$Nkd[,min(mdist)])-0.01
  registerDoParallel(cores=ncpus)
  ops.count = foreach (i=1:squares[,.N], .packages=c("data.table","rstan"), .combine="rbind") %dopar% {
    #extract relevant portion of data
    extracted = get_cs_subset(counts, biases.aug, begin1=squares[i,begin1], end1=squares[i,end1],
                              begin2=squares[i,begin2], end2=squares[i,end2], fill.zeros=T)
    #run fit
    a=system.time(op <- optimize_extradiag(smfitex, extracted$biases1, extracted$biases2, extracted$counts, dmin, dmax,
                                           bf_per_decade = bf_per_decade, verbose = verbose, iter = iter))
    #report data
    setkey(extracted$biases1, id, pos)
    setkey(extracted$counts, id1, pos1)
    extracted$counts[,c("log_nu1","log_delta1"):=extracted$biases1[extracted$counts,.(log_nu,log_delta)]]
    nu1off=mean(extracted$counts[,log_nu1])
    delta1off=mean(extracted$counts[,log_delta1])
    setkey(extracted$biases2, id, pos)
    setkey(extracted$counts, id2, pos2)
    extracted$counts[,c("log_nu2","log_delta2"):=extracted$biases2[extracted$counts,.(log_nu,log_delta)]]
    nu2off=mean(extracted$counts[,log_nu2])
    delta2off=mean(extracted$counts[,log_delta2])
    center=extracted$counts[,min(distance)+max(distance)]/2
    border=extracted$counts[,max(distance)-min(distance)]/2+1
    #
    setkey(extracted$counts, id1,id2,pos1,pos2,distance)
    ret=extracted$counts[,.(id1=id1+extracted$beginrange1-1,id2=id2+extracted$beginrange2-1,
                            pos1,pos2,distance,weight=(border-abs(distance-center))**2)]
    ret[,c("square.begin1","square.begin2","log_decay_close","log_decay_far","log_decay_up","log_decay_down"):=list(
           squares[i,begin1], squares[i,begin2],
           op$par$eC+op$par$log_decay_close+nu1off-delta1off+nu2off+delta2off,
           op$par$eC+op$par$log_decay_far  +nu1off+delta1off+nu2off-delta2off,
           op$par$eC+op$par$log_decay_up   +nu1off+delta1off+nu2off+delta2off,
           op$par$eC+op$par$log_decay_down +nu1off-delta1off+nu2off-delta2off)]
    ret
  }
  ### reconstruct count estimates: log_decay and beta_diag_centered
  message("*** reconstruct diagonal decay")
  #qplot(ops.count[,id],binwidth=1)
  #ops.count[,cat:=paste(square.begin1,square.begin2)]
  #ops.count[,dbin:=factor(square.begin2-square.begin1)]
  #ggplot(ops.count)+geom_line(aes(distance,weight,colour=cat))+guides(colour=F)+facet_grid(dbin~.)
  #ggplot(ops.count[weight>1e9])+geom_line(aes(distance,log_decay_far,colour=cat))+guides(colour=F)+facet_grid(.~dbin)
  #ggplot(ops.count[square.begin1==square.begin2&square.begin1==3189.0])+geom_line(aes(distance,log_decay_up))
  #
  stepsz=1/distance_bins_per_decade
  #if (circularize>0) {
  #  dmax=biases.aug[,max(pmin(pos,circularize-pos+1))]
  #  dmin=biases.aug[,min(pmin(pos,circularize-pos+1))]
  #  dbins=10**seq(0,log10(dmax-dmin)+stepsz,stepsz)
  #} else {
    dbins=10**seq(0,biases.aug[,log10(max(pos)-min(pos))]+stepsz,stepsz)
  #}
  ops.count[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #ggplot(ops.count)+geom_line(aes(distance,log_decay_far,group=cat,colour=bdist))+guides(colour=F)+facet_grid(.~dbin)
  test=ops.count[,.(log_decay_close=weighted.mean(log_decay_close, weight),
             log_decay_far=weighted.mean(log_decay_far, weight),
             log_decay_up=weighted.mean(log_decay_up, weight),
             log_decay_down=weighted.mean(log_decay_down, weight),
             ncounts=.N), keyby=bdist] #could only take one decay, and even weights since eC is estimated later
  bin.begin=test[,bdist]
  bin.end=test[,bdist]
  levels(bin.begin) <- tstrsplit(as.character(levels(bin.begin)), "[[,]")[2][[1]]
  levels(bin.end) <- tstrsplit(as.character(levels(bin.end)), "[[,)]")[2][[1]]
  test[,begin:=as.integer(as.character(bin.begin))]
  test[,end:=as.integer(as.character(bin.end))]
  smfitdec=stan_model("cs_norm_fit_decay_params.stan")
  op=optimize_decay_params(smfitdec, test, counts, dmin, dmax, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)
  test[,log_decay:=op$par$log_decay+op$par$eC]
  retlist$beta_diag_centered=op$par$beta_diag_centered
  ### reconstruct count estimates: eC
  message("*** reconstruct count exposure")
  smfiteC=stan_model("cs_norm_fit_eC.stan")
  op=optimize_eC(smfiteC, biases, counts[sample(.N,min(.N,100000))], retlist, dmin, dmax, bf_per_decade=bf_per_decade, verbose=verbose)
  retlist$eC=op$par$eC
  retlist$alpha=op$par$alpha
  #ggplot(melt(test, id.vars=c("bdist","begin","end")))+geom_point(aes(bdist,value,colour=variable))
  return(list(par=retlist))
}

postprocess = function(biases, counts, meanfield, op, resolution=10000, ncores=30, predict.all.means=T) {
  smpred = stan_model(file = "cs_norm_predict.stan")
  smbin = stan_model("cs_norm_predict_binned.stan")
  smdisp = stan_model("cs_norm_binned_dispersions.stan")
  ### run remaining steps
  if (predict.all.means==T) {
    message("*** predict all means")
    op$pred=predict_all_meanfield(smpred, biases, counts, meanfield, op, verbose=T)$par
  }
  message("*** buid binned matrices")
  op$binned=get_binned_matrices(smbin, biases, counts, meanfield, op, resolution=resolution)
  message("*** estimate dispersions")
  op$disp=get_dispersions(smdisp, op$binned$mat)$par
  message("*** detect interactions")
  op$mat=detect_interactions(op$binned$mat, op$disp$dispersion, ncores=ncores) #interaction detection using binned dispersion estimates
  return(op)
}


#read_and_prepare("/scratch/rao/mapped/HICall_both_filled_map_chr1.tsv", "data/rao_HICall_chr1_all", skip="SRR")

#read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
#                 "data/caulo_BglIIr1_all", skip="SRR", circularize=4042929)
#read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate2_reads_int.tsv",
#                 "data/caulo_BglIIr2_all", skip="SRR", circularize=4042929)
#read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
#                 "data/caulo_NcoI_all", skip="SRR", circularize=4042929)

#read_and_prepare("/scratch/ralph/HiC/3_Mapped/Bcell_Sox2_10Mb_both_filled_map.tsv", "data/ralph_Bcell_Sox2", skip="HWI")
#read_and_prepare("/scratch/ralph/HiC/3_Mapped/EScell_Sox2_10Mb_both_filled_map.tsv", "data/ralph_EScell_Sox2", skip="HWI")

dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = 0, eRJ = 5, eDE = 7)
counts=dset$counts
biases=dset$biases
dset_statistics(biases,counts)
meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)

prefix="caulo_NcoI_all"
#prefix="ralph_Bcell_Sox2"
biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))
#meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)
#save(meanfield, file = paste0("data/",prefix,"_meanfield_100.RData"))
load(paste0("data/",prefix,"_meanfield_100.RData"))

#ralph: restrict to 30500000,32500000
biases=biases[pos>=30500000&pos<=30600000]
#biases=biases[pos>=30500000&pos<=32500000]
beginrange=biases[1,id]
endrange=biases[.N,id]
biases[,id:=id-beginrange+1]
prefix=biases[,paste0(prefix,"_",min(pos),"-",max(pos))]
counts=counts[id1>=beginrange&id1<=endrange&id2>=beginrange&id2<=endrange]
counts[,c("id1","id2"):=list(id1-beginrange+1,id2-beginrange+1)]
meanfield$Nkl=meanfield$Nkl[id>=beginrange&id<=endrange][,.(id=id-beginrange+1,count,N)]
meanfield$Nkr=meanfield$Nkr[id>=beginrange&id<=endrange][,.(id=id-beginrange+1,count,N)]

#caulo/toy: restrict to 100 consecutive rsites
beginrange=1
endrange=975
beginrange=1
endrange=200
beginrange=1
endrange=80
beginrange=22
endrange=96
beginrange=49
endrange=119
biases=biases[beginrange:endrange]
biases[,id:=id-beginrange+1]
prefix=biases[,paste0("caulo_NcoI_",min(pos),"-",max(pos))]
counts=counts[id1>=beginrange&id1<=endrange&id2>=beginrange&id2<=endrange]
counts[,c("id1","id2"):=list(id1-beginrange+1,id2-beginrange+1)]
meanfield$Nkl=meanfield$Nkl[id>=beginrange&id<=endrange][,.(id=id-beginrange+1,count,N)]
meanfield$Nkr=meanfield$Nkr[id>=beginrange&id<=endrange][,.(id=id-beginrange+1,count,N)]


counts=fill_zeros(counts,biases)
counts[,distance:=abs(pos2-pos1)] #not filled
write.table(biases, file = paste0(prefix,"_biases.dat"), quote=F, row.names = F)
write.table(counts, file = paste0(prefix,"_counts.dat"), quote=F, row.names = F)
save(meanfield, file = paste0(prefix, "_meanfield_100.RData"))



#### optimization wihout prior guesses
smfit = stan_model(file = "cs_norm_fit.stan")
smpred = stan_model(file = "cs_norm_predict.stan")
smbin = stan_model("cs_norm_predict_binned.stan")
smdisp = stan_model("cs_norm_binned_dispersions.stan")
maxcount=-1
a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                           bf_per_decade=5, verbose = T, iter=100000)) #, tol_rel_grad=1e3, tol_rel_obj=1e3
op$pred=predict_all_meanfield(smpred, biases, counts, meanfield, op, verbose=T)$par
op$binned=get_binned_matrices(smbin, biases, counts, meanfield, op, resolution=10000)
op$disp=get_dispersions(smdisp, op$binned$mat)$par
op$mat=detect_interactions(op$binned$mat, op$disp$dispersion, ncores=30) #interaction detection using binned dispersion estimates
mat=op$mat
mat2=detect_interactions(op$binned$mat, op$par$alpha*op$binned$mat[,ncounts], ncores=1) #interaction detection using original dispersion estimate
save(op, file=paste0("data/",prefix,"_op_maxcount_",maxcount,".RData"))


#plot observed / expected matrices and histograms
ggplot(data.table(dist=op$binned$dist, decay=op$binned$decay))+geom_point(aes(dist,decay))+scale_x_log10()+scale_y_log10()
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(ncounts)))
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(observed)))
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(expected)))
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=mat[is.interaction==T])
ggplot(mat2)+geom_raster(aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=mat[is.interaction==T])
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=lFC))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=mat[is.interaction==T])
#ggplot(mat)+geom_raster(aes(begin1,begin2,fill=prob.observed.gt.expected>0.95))
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=is.interaction))
ggplot(mat2)+geom_raster(aes(begin1,begin2,fill=is.interaction))
ggplot(mat)+geom_histogram(aes(lFC))


#compare dispersion estimates
message("estimation during fit ",op$par$alpha*op$binned$mat[,mean(ncounts)])
message("estimation on binned matrices ",mean(op$disp$dispersion))
message("estimation with glm.nb ", op$binned$mat[,theta.ml(observed,expected)])
ggplot(data.table(dfit=op$par$alpha*op$binned$mat[,ncounts],dbin=op$disp$dispersion))+geom_point(aes(dbin,dfit))+stat_function(fun=identity)

#TODO:
#- estimate dispersion on each bin individually
#- estimate dispersion in the form alpha_i * alpha_j
#- implement peak calling

ggplot(data.table(x=c(132,133)),aes(x))+
  stat_function(fun=function(x){dgamma(x,7089163,rate=53574)},colour="blue")+
  stat_function(fun=function(x){dgamma(x,7089277,rate=53575)},colour="red")

a=thresholds_estimator(114,132,100000)

#### effect of fitting parts of a matrix
load("data/caulo_NcoI_3189-363550_op_maxcount_-1.RData")
opall=op
load("data/caulo_NcoI_3189-147895_op_maxcount_-1.RData")
opfirst=op
load("data/caulo_NcoI_52281-194220_op_maxcount_-1.RData")
opsecond=op
load("data/caulo_NcoI_100417-249556_op_maxcount_-1.RData")
opthird=op

ggplot()+geom_raster(data=opfirst$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opfirst$mat[is.interaction==T])+
  geom_raster(data=opsecond$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opsecond$mat[is.interaction==T])
#
ggplot(opall$mat)+geom_raster(aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
#
ggplot()+geom_raster(data=opfirst$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opfirst$mat[is.interaction==T])+
  geom_raster(data=opall$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
#
ggplot()+geom_raster(data=opsecond$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opsecond$mat[is.interaction==T])+
  geom_raster(data=opall$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
#
ggplot()+geom_raster(data=opthird$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opthird$mat[is.interaction==T])+
  geom_raster(data=opall$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
ggsave("images/caulo_NcoI_3189-363550_op_all_vs_third.png", width=10, height=7.5)
#
biases[1:200,log_nu_all:=opall$par$log_nu]
biases[1:80,log_nu_first:=opfirst$par$log_nu]
biases[22:96,log_nu_second:=opsecond$par$log_nu]
ggplot(biases[1:100])+geom_point(aes(pos,log_nu_all,colour="all"))+geom_point(aes(pos,log_nu_first,colour="first"))+geom_point(aes(pos,log_nu_second,colour="second"))
ggplot(biases[1:100])+geom_point(aes(pos,abs(log_nu_first-log_nu_all),colour="first"))+geom_point(aes(pos,abs(log_nu_second-log_nu_all),colour="second"))
ggplot(biases[1:100])+geom_point(aes(pos,abs(1-exp(log_nu_first-log_nu_all)),colour="first"))+geom_point(aes(pos,abs(1-exp(log_nu_second-log_nu_all)),colour="second"))
biases[1:200,log_delta_all:=opall$par$log_delta]
biases[1:80,log_delta_first:=opfirst$par$log_delta]
biases[22:96,log_delta_second:=opsecond$par$log_delta]
ggplot(biases[1:100])+geom_point(aes(pos,log_delta_all,colour="all"))+geom_point(aes(pos,log_delta_first,colour="first"))+geom_point(aes(pos,log_delta_second,colour="second"))
ggplot(biases[1:100])+geom_point(aes(pos,abs(log_delta_first-log_delta_all),colour="first"))+geom_point(aes(pos,abs(log_delta_second-log_delta_all),colour="second"))
ggplot(biases[1:100])+geom_point(aes(pos,abs(1-exp(log_delta_first-log_delta_all)),colour="first"))+geom_point(aes(pos,abs(1-exp(log_delta_second-log_delta_all)),colour="second"))
ggplot()+scale_x_log10()+scale_y_log10()+
  geom_line(data=data.table(x=opall$binned$distance,y=opall$binned$decay),aes(x,y,colour="all"))+
  geom_line(data=data.table(x=opfirst$binned$distance,y=opfirst$binned$decay),aes(x,y,colour="first"))+
  geom_line(data=data.table(x=opsecond$binned$distance,y=opsecond$binned$decay),aes(x,y,colour="second"))

#### effect of mean field
load("data/caulo_NcoI_3189-363550_op_maxcount_-1.RData")
opall=op
load("data/caulo_NcoI_3189-363550_op_maxcount_0.RData")
opmf0=op
load("data/caulo_NcoI_3189-363550_op_maxcount_1.RData")
opmf1=op
load("data/caulo_NcoI_3189-363550_op_maxcount_4.RData")
opmf4=op

ggplot()+geom_raster(data=opmf0$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opmf0$mat[is.interaction==T])+
  geom_raster(data=opall$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
ggplot()+geom_raster(data=opmf1$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opmf1$mat[is.interaction==T])+
  geom_raster(data=opall$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
ggplot()+geom_raster(data=opmf4$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=opmf4$mat[is.interaction==T])+
  geom_raster(data=opall$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opall$mat[is.interaction==T])
ggsave("images/caulo_NcoI_3189-363550_op_maxcount_-1_vs_4.png", width=10, height=7.5)


### effect of parallelization
#opall <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=-1, bf_per_kb=1,
#                                bf_per_decade=5, verbose = T, iter=100000)
#opall=postprocess(biases, counts, meanfield, opall, resolution=10000, ncores=30)


coverage=4
square.size=100000
oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=1,
                      bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=100000, ncpus=30, homogenize=F)
oppar=postprocess(biases, counts, meanfield, oppar, resolution=10000, ncores=30)
save(oppar, file = paste0("data/caulo_NcoI_3189-363550_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k.RData"))

load("data/caulo_NcoI_3189-363550_op_maxcount_-1.RData")
opserial=op
load("data/caulo_NcoI_3189-363550_op_maxcount_-1_parallel_inhomogeneous_cov10X_sq150k.RData")

#lFC histogram and matrices
ggplot(rbind(opserial$mat[,.(bin1,bin2,dset="serial",lFC)],oppar$mat[,.(bin1,bin2,dset="parallel",lFC)]))+geom_histogram(aes(lFC,fill=dset),position="dodge")
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_hist.png"), width=10, height=7.5)
ggplot()+geom_raster(data=oppar$mat,aes(begin1,begin2,fill=lFC))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+
  geom_raster(data=opserial$mat,aes(begin2,begin1,fill=lFC))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opserial$mat[is.interaction==T])
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_mat.png"), width=10, height=7.5)
#normalized matrix
ggplot()+geom_raster(data=oppar$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+
  geom_raster(data=opserial$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opserial$mat[is.interaction==T])
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_normalized.png"), width=10, height=7.5)
#genomic biases
ggplot(melt(data.table(id=biases[,pos],serial=opserial$par$log_nu,parallel=oppar$par$log_nu),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu.png"), width=10, height=7.5)
ggplot(melt(data.table(id=biases[,pos],serial=opserial$par$log_delta,parallel=oppar$par$log_delta),id.vars="id",variable.name="delta"))+geom_line(aes(id,value,colour=delta))
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_delta.png"), width=10, height=7.5)
ggplot(rbind(data.table(bias="nu",serial=opserial$par$log_nu,parallel=oppar$par$log_nu),data.table(bias="delta",serial=opserial$par$log_delta,parallel=oppar$par$log_delta)))+
  geom_point(aes(serial,parallel))+facet_grid(.~bias)+stat_function(fun=identity)
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu_delta.png"), width=10, height=7.5)
#diagonal biases
ggplot(melt(data.table(dist=counts[,distance],serial=opserial$pred$log_decay_close,par=oppar$pred$log_decay_far, key="dist"),id.vars="dist"))+
  geom_line(aes(dist,value,colour=variable))+scale_x_log10()
ggsave(filename=paste0("images/caulo_NcoI_3189-363550_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_decay.png"), width=10, height=7.5)
c(opserial$mat[,mean(lFC)], oppar$mat[,mean(lFC)])
c(opserial$par$eRJ,oppar$par$eRJ)
c(opserial$par$eDE,oppar$par$eDE)
c(opserial$par$eC,oppar$par$eC)
c(opserial$disp$alpha,oppar$disp$alpha)




### effect of parallelization
coverage=4
square.size=150000
oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=1,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=100000, ncpus=30, homogenize=F)
oppar=postprocess(biases, counts, meanfield, oppar, resolution=30000, ncores=30, predict.all.means=F)
oppar$ice=iterative_normalization(oppar$mat, niterations=1, resolution=30000, return.binned=T)
save(oppar, file = paste0("data/caulo_NcoI_",biases[,min(pos)],"-",biases[,max(pos)],"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k.RData"))

setkey(oppar$mat, begin1,begin2)
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(normalized)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none") 
ggsave(filename = paste0("images/caulo_NcoI_",biases[,min(pos)],"-",biases[,max(pos)],"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_normalized.png"), width=10, height=7.5)
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=lFC))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=lFC))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+theme(legend.position = "none") 
ggsave(filename = paste0("images/caulo_NcoI_",biases[,min(pos)],"-",biases[,max(pos)],"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_mat.png"), width=10, height=7.5)
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(observed)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+scale_fill_gradient(low="white", high="black")+theme(legend.position = "none") 
ggsave(filename = paste0("images/caulo_NcoI_",biases[,min(pos)],"-",biases[,max(pos)],"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_observed.png"), width=10, height=7.5)
ggplot()+geom_raster(data=oppar$ice, aes(begin1,begin2,fill=log(N)))+geom_raster(data=oppar$ice, aes(begin2,begin1,fill=log(N)))+scale_fill_gradient(low="white", high="black")+  theme(legend.position = "none") 
ggsave(filename = paste0("images/caulo_NcoI_",biases[,min(pos)],"-",biases[,max(pos)],"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_ICE.png"), width=10, height=7.5)+  theme(legend.position = "none") 
