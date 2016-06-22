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

bin_counts = function(counts, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL) {
  if (is.null(b1)) b1=counts[,min(pos1)]-1
  if (is.null(b2)) b2=counts[,min(pos2)]-1
  if (is.null(e1)) e1=counts[,max(pos1)]-1
  if (is.null(e2)) e2=counts[,max(pos2)]-1
  bins1=seq(b1,counts[,max(pos1)]+resolution,resolution)
  bins2=seq(b2,counts[,max(pos2)]+resolution,resolution)
  mcounts=melt(counts,measure.vars=c("contact.close","contact.far","contact.up","contact.down"),
               variable.name = "category", value.name = "count")[count>0]
  #
  sub = mcounts[,.(pos1,pos2,bin1=cut2(pos1, bins1, oneval=F, onlycuts=T, digits=10),
                   bin2=cut2(pos2, bins2, oneval=F, onlycuts=T, digits=10), category, count)
                ][,.(N=sum(count)),by=c("bin1","bin2")]
  #
  sub[,begin1:=do.call(as.integer, tstrsplit(as.character(bin1), "[[,]")[2])]
  sub[,end1:=do.call(as.integer, tstrsplit(as.character(bin1), "[],)]")[2])]
  sub[,begin2:=do.call(as.integer, tstrsplit(as.character(bin2), "[[,]")[2])]
  sub[,end2:=do.call(as.integer, tstrsplit(as.character(bin2), "[],)]")[2])]
  #
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
  enrich=enrich[,.(N=sum(N)),by=c("re.closest1","re.closest2","category")]
  #
  #get rsites that a) have counts on them and b) have at least one dangling/rejoined end.
  #message("RSites list")
  #rsites.count=enrich[category %in% c("contact.close", "contact.down", "contact.far", "contact.up"),
  #                    unique(c(re.closest1,re.closest2))]
  #rsites.bias=enrich[category %in% c("dangling.L", "dangling.R", "rejoined"),
  #                   unique(c(re.closest1,re.closest2))]
  #rsites=data.table(re.pos=intersect(rsites.count, rsites.bias))
  rsites = data.table(re.pos=enrich[,unique(c(re.closest1,re.closest2))])
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
  X=rbind(X,rsites[id%nin%X[,id],.(id,pos=re.pos,dangling.L=0,dangling.R=0,rejoined=0)])
  setkey(X,id)
  stopifnot(all(X[,id==.I])) #IDs must start at one and have no gaps, for current stan implementation
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

csnorm_fit = function(model, biases, counts, dmin, dmax, bf_per_kb=1, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  dmax=counts[,max(distance)]+0.01
  dmin=counts[,min(distance)]-0.01
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far], counts_up=counts[,contact.up], counts_down=counts[,contact.down])
  if (verbose==T) {
    message("CS norm: fit")
    message("Krow        : ", Krow)
    message("Kdiag       : ", Kdiag)
    message("Biases      : ", biases[,.N])
    message("Counts      : ", 4*counts[,.N])
    message("% zeros     : ", (counts[contact.close==0,.N]+counts[contact.far==0,.N]+counts[contact.up==0,.N]+counts[contact.down==0,.N])/counts[,.N]*100)
  }
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

csnorm_fit_extradiag = function(model, biases1, biases2, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( S1=biases1[,.N], S2=biases2[,.N], Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=counts[,.N], cidx=t(data.matrix(counts[,.(id1,id2)])), dist=counts[,distance],
               counts_close=counts[,contact.close], counts_far=counts[,contact.far], counts_up=counts[,contact.up], counts_down=counts[,contact.down],
               log_nu1=biases1[,log_nu], log_nu2=biases2[,log_nu], log_delta1=biases1[,log_delta], log_delta2=biases2[,log_delta])
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

get_cs_subset = function(counts, biases, begin1, end1, begin2=NULL, end2=NULL, fill.zeros=T, circularize=-1) {
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
    if (circularize>0) {
      counts.local[,distance:=pmin(abs(pos2-pos1), circularize+1-abs(pos2-pos1))]
    } else {
      counts.local[,distance:=abs(pos2-pos1)]
    }
  }
  return(list(counts=counts.local,biases1=biases.1,biases2=biases.2,
              beginrange1=beginrange1, endrange1=endrange1, beginrange2=beginrange2, endrange2=endrange2))
}

run_split_parallel_squares = function(biases, square.size, coverage, diag.only=F) {
  minpos=biases[,min(pos)]
  maxpos=biases[,max(pos)+1]
  true.square.size=(maxpos-minpos)/round((maxpos-minpos)/square.size)
  begins=seq(minpos,maxpos-true.square.size,by=true.square.size/coverage)
  if (diag.only==T) {
    squares=SJ(begins,begins)
  } else {
    squares=CJ(begins,begins)[V2>=V1]
  }
  setnames(squares, c("begin1","begin2"))
  squares[,c("end1","end2"):=list(begin1+true.square.size,begin2+true.square.size)]
  diagsquares=squares[begin1==begin2]
  #ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2),alpha=0.1)+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=factor((begin2-begin1)/true.square.size)),alpha=0.1)+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin2-end1,xmax=end2-begin1,ymin=sqid,ymax=sqid+1),alpha=0.1)+scale_x_log10()
  return(list(squares=squares,diagsquares=diagsquares,true.square.size=true.square.size))
}

run_split_parallel_biases = function(smfit, counts, biases, begin, end, dmin, dmax, bf_per_kb, bf_per_decade, verbose, iter, outprefix=NULL, circularize=-1) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases, begin1=begin, end1=end, fill.zeros=T, circularize=circularize)
  #run fit
  a=system.time(output <- capture.output(op <- csnorm_fit(smfit, extracted$biases1, extracted$counts, dmin, dmax,
                                                          bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)))
  op$runtime=a[1]+a[4]
  op$output=output
  if (!is.null(outprefix)) {
    save(op, file=paste0(outprefix, "_biases_op_",begin,".RData"))
  }
  #report data
  center=(end-begin)/2
  ret=extracted$biases1[,.(id=id+extracted$beginrange1-1,pos,weight=(center-abs(pos-center))**2)]
  ret[,c("log_mean_DL", "log_mean_DR", "log_mean_RJ"):=list(op$par$log_mean_DL, op$par$log_mean_DR, op$par$log_mean_RJ)]
  ret[,square.begin:=begin]
  ret[,c("lambda_delta","lambda_nu","nsites","alpha"):=list(op$par$lambda_delta,op$par$lambda_nu,extracted$biases1[,.N],op$par$alpha)]
  if (!is.null(outprefix)) {
    save(ret, file=paste0(outprefix, "_biases_ret_",begin,".RData"))
  }
  return(list(ret=ret, out=tail(op$output,1), runtime=op$runtime))
}

run_split_parallel_biases_homogenize = function(model, dt, lambda_nu, lambda_delta, bf_per_kb=1, iter=10000, verbose=T, init=0, ...) {
  Krow=round(dt[,(max(pos)-min(pos))/1000*bf_per_kb])
  data = list( Krow=Krow, S=dt[,.N], cutsites=dt[,pos],
               log_mean_RJ=dt[,log_mean_RJ], log_mean_DL=dt[,log_mean_DL], log_mean_DR=dt[,log_mean_DR],
               lambda_nu=lambda_nu, lambda_delta=lambda_delta)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

run_split_parallel_counts = function(model, counts, biases.aug, begin1, end1, begin2, end2, dmin, dmax,
                                     bf_per_decade=5, verbose=T, iter=100000, outprefix=NULL, circularize=-1L) {
  #extract relevant portion of data
  extracted = get_cs_subset(counts, biases.aug, begin1=begin1, end1=end1,
                            begin2=begin2, end2=end2, fill.zeros=T, circularize=circularize)
  #run fit
  a=system.time(output <- capture.output(op <- csnorm_fit_extradiag(model, extracted$biases1, extracted$biases2, extracted$counts,
                                                                    dmin, dmax, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)))
  op$runtime=a[1]+a[4]
  op$output=output
  if (!is.null(outprefix)) {
    save(op, file=paste0(outprefix, "_counts_op_",begin1,"_",begin2,".RData"))
  }
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
    begin1, begin2,
    op$par$eC+op$par$log_decay+nu1off-delta1off+nu2off+delta2off,
    op$par$eC+op$par$log_decay+nu1off+delta1off+nu2off-delta2off,
    op$par$eC+op$par$log_decay+nu1off+delta1off+nu2off+delta2off,
    op$par$eC+op$par$log_decay+nu1off-delta1off+nu2off-delta2off)]
  if (!is.null(outprefix)) {
    save(ret, file=paste0(outprefix, "_counts_ret_",begin1,"_",begin2,".RData"))
  }
  return(list(ret=ret, out=tail(op$output,1), runtime=op$runtime))
}

run_split_parallel_counts_decay = function(model, dt, counts, dmin, dmax, bf_per_decade=5, iter=10000, verbose=T, init=0, ...) {
  dist=dt[,(begin+end)/2]
  fij=dt[,(log_decay_close+log_decay_far+log_decay_up+log_decay_down)/4]
  ncounts=dt[,ncounts]
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, dmin=dmin, dmax=dmax,
               N=length(dist), dist=dist, fij=fij, ncounts=ncounts)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...)
}

run_split_parallel_counts_eC = function(model, biases, counts, retlist, dmin, dmax, bf_per_decade=5, verbose=T, iter=1000) {
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round((log10(dmax)-log10(dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=dmin, dmax=dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               log_nu=retlist$log_nu, log_delta=retlist$log_delta, beta_diag_centered=retlist$beta_diag_centered)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)
}

run_split_parallel = function(counts, biases, square.size=100000, coverage=4, coverage.extradiag=1, bf_per_kb=1, bf_per_decade=5,
                              distance_bins_per_decade=100, verbose = F, iter=100000, ncpus=30, homogenize=F, outprefix=NULL, circularize=-1L,
                              ops.bias=NULL, ops.count=NULL) {
  stopifnot(!(circularize==-1 && counts[,max(distance)]<biases[,max(pos)-min(pos)]))
  ### build squares
  message("*** build squares")
  retval.squares = run_split_parallel_squares(biases, square.size, coverage, diag.only=T)
  diagsquares=retval.squares$diagsquares
  true.square.size=retval.squares$true.square.size
  retval.squares = run_split_parallel_squares(biases, square.size, coverage.extradiag, diag.only=F)
  stopifnot(true.square.size==retval.squares$true.square.size)
  squares=retval.squares$squares
  message("Run in parallel")
  message("   square size: ",true.square.size)
  message("   diagonal coverage: ",coverage)
  message("   extradiagonal coverage: ",coverage.extradiag)
  message("   ",squares[,.N], " squares")
  message("   ",diagsquares[,.N], " on the diagonal")
  ### fit genomic biases using squares close to diagonal
  message("*** fit genomic biases")
  dmax=counts[,max(distance)]+0.01
  dmin=counts[,min(distance)]-0.01
  registerDoParallel(cores=ncpus)
  smfit = stan_model(file = "cs_norm_fit.stan")
  output.binder = function(x) {
    a=sapply(names(x[[1]]), function(i) sapply(x, "[[", i,simplify=F))
    list(par=rbindlist(a[,1]), out=as.character(a[,2]), runtime=as.numeric(a[,3]))
  }
  if (is.null(ops.bias)) {
    ops.bias = foreach (begin=diagsquares[,begin1], end=diagsquares[,end1], .packages=c("data.table","rstan")) %dopar% 
      run_split_parallel_biases(smfit, counts, biases, begin, end, dmin, dmax, bf_per_kb = bf_per_kb,
                                bf_per_decade = bf_per_decade, verbose = verbose, iter = iter, outprefix=outprefix, circularize=circularize)
    ops.bias = output.binder(ops.bias)
    if (!is.null(outprefix)) {save(ops.bias, file=paste0(outprefix, "_ops_bias.RData"))}
  }
  ### reconstruct bias estimates: eRJ eDE log_nu and log_delta
  message("*** reconstruct genomic biases")
  info=ops.bias$par[,.SD[1],by=square.begin,.SDcols=c("lambda_nu","lambda_delta","nsites","alpha")]
  #ggplot(info)+geom_line(aes(square.begin,lambda_nu))+scale_y_log10()+geom_hline(yintercept=info[,exp(median(log(lambda_nu)))])
  means=ops.bias$par[,.(log_mean_RJ=weighted.mean(log_mean_RJ, weight),
                        log_mean_DL=weighted.mean(log_mean_DL, weight),
                        log_mean_DR=weighted.mean(log_mean_DR, weight),
                        ncounts=.N), keyby=c("id","pos")]
  if (homogenize==T) {
    message("*** homogenize genomic biases")
    smfitgen=stan_model("cs_norm_fit_genomic_params.stan")
    op=run_split_parallel_biases_homogenize(smfitgen, means, lambda_nu=info[,exp(mean(log(lambda_nu)))], lambda_delta=info[,exp(mean(log(lambda_delta)))],
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
    retlist=list(eRJ=ops.bias$par[,mean(log_mean_RJ)], eDE=ops.bias$par[,mean(log_mean_DL+log_mean_DR)/2])
    biases.aug=means[biases]
    stopifnot(!any(is.na(biases.aug)))
    biases.aug[,log_nu:=log_mean_RJ-retlist$eRJ]
    biases.aug[,log_delta:=(log_mean_DL-log_mean_DR)/2]
    retlist$log_nu=biases.aug[,log_nu]
    retlist$log_delta=biases.aug[,log_delta]
  }
  retlist$out.bias=ops.bias$out
  retlist$runtime.bias=ops.bias$runtime
  #ggplot()+geom_line(data=data.table(x=biases.aug[,pos],y=retlist$log_nu),aes(x,y))
  #
  ### fit remaining data
  message("*** fit diagonal decay")
  smfitex=stan_model("cs_norm_fit_extradiag.stan")
  registerDoParallel(cores=ncpus)
  if (is.null(ops.count)) {
    ops.count = foreach (begin1=squares[,begin1], end1=squares[,end1], begin2=squares[,begin2], end2=squares[,end2],
                         .packages=c("data.table","rstan")) %dopar% 
      run_split_parallel_counts(smfitex, counts, biases.aug, begin1, end1, begin2, end2, dmin, dmax,
                                bf_per_decade = bf_per_decade, verbose = verbose, iter = iter, outprefix=outprefix, circularize=circularize)
    ops.count=output.binder(ops.count)
    if (!is.null(outprefix)) {save(ops.count, file=paste0(outprefix, "_ops_count.RData"))}
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
  dbins=10**seq(log10(dmin),log10(dmax)+stepsz,stepsz)
  #}
  ops.count$par[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=5)]
  #ops.count$par[,grid:=factor(square.begin2-square.begin1)]
  #ggplot(ops.count$par[sample(.N,100000)])+geom_line(aes(distance,log_decay_far,group=square.begin1-square.begin2,colour=bdist))+guides(colour=F)+facet_grid(.~grid)
  test=ops.count$par[,.(log_decay_close=weighted.mean(log_decay_close, weight),
                        log_decay_far=weighted.mean(log_decay_far, weight),
                        log_decay_up=weighted.mean(log_decay_up, weight),
                        log_decay_down=weighted.mean(log_decay_down, weight),
                        ncounts=.N), keyby=bdist] #could only take one decay, and even weights since eC is estimated later
  bin.begin=test[,bdist]
  bin.end=test[,bdist]
  levels(bin.begin) <- tstrsplit(as.character(levels(bin.begin)), "[[,]")[2][[1]]
  levels(bin.end) <- tstrsplit(as.character(levels(bin.end)), "[],)]")[2][[1]]
  test[,begin:=as.integer(as.character(bin.begin))]
  test[,end:=as.integer(as.character(bin.end))]
  smfitdec=stan_model("cs_norm_fit_decay_params.stan")
  op=run_split_parallel_counts_decay(smfitdec, test, counts, dmin, dmax, bf_per_decade = bf_per_decade, verbose = verbose, iter = iter)
  test[,log_decay:=op$par$log_decay+op$par$eC]
  retlist$beta_diag_centered=op$par$beta_diag_centered
  ### reconstruct count estimates: eC
  message("*** reconstruct count exposure")
  smfiteC=stan_model("cs_norm_fit_eC.stan")
  op=run_split_parallel_counts_eC(smfiteC, biases, counts[sample(.N,min(.N,100000))], retlist, dmin, dmax, bf_per_decade=bf_per_decade, verbose=verbose)
  retlist$eC=op$par$eC
  retlist$alpha=op$par$alpha
  #ggplot(test)+geom_point(aes(bdist, log_decay))
  #ggplot(melt(test, id.vars=c("bdist","begin","end")))+geom_point(aes(bdist,value,colour=variable))
  return(list(par=retlist, out.bias=ops.bias$out, out.count=ops.count$out, runtime.count=ops.count$runtime, runtime.bias=ops.bias$runtime, dmin=dmin, dmax=dmax))
}

run_split_parallel_recovery = function(counts, biases, outprefix, square.size=100000, coverage=4, coverage.extradiag=1, bf_per_kb=1, bf_per_decade=5,
                                       distance_bins_per_decade=100, verbose = F, iter=100000, ncpus=30, homogenize=F, circularize=-1L) {
  output.binder = function(x) {
    a=sapply(names(x[[1]]), function(i) sapply(x, "[[", i,simplify=F))
    list(par=rbindlist(a[,1]), out=as.character(a[,2]), runtime=as.numeric(a[,3]))
  }
  ops.bias=sapply(X=Sys.glob(paste0(outprefix,"_biases_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
  ops.bias = output.binder(ops.bias)
  ops.count=sapply(X=Sys.glob(paste0(outprefix,"_counts_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
  ops.count = output.binder(ops.count)
  run_split_parallel(counts, biases, square.size, coverage, coverage.extradiag, bf_per_kb, bf_per_decade,
                     distance_bins_per_decade, verbose, iter, ncpus, homogenize, outprefix, circularize,
                     ops.bias, ops.count)
}

csnorm_predict_all = function(model, biases, counts, opt, bf_per_decade=5, verbose=T) {
  cclose=counts[,.(id1,id2,distance,count=contact.close)]
  cfar=counts[,.(id1,id2,distance,count=contact.far)]
  cup=counts[,.(id1,id2,distance,count=contact.up)]
  cdown=counts[,.(id1,id2,distance,count=contact.down)]
  Kdiag=round((log10(opt$dmax)-log10(opt$dmin))*bf_per_decade)
  data = list( Kdiag=Kdiag, S=biases[,.N], cutsites=biases[,pos], dmin=opt$dmin, dmax=opt$dmax,
               Nclose=cclose[,.N], counts_close=cclose[,count], index_close=t(data.matrix(cclose[,.(id1,id2)])), dist_close=cclose[,distance],
               Nfar=cfar[,.N],     counts_far=cfar[,count],     index_far=t(data.matrix(cfar[,.(id1,id2)])), dist_far=cfar[,distance],
               Nup=cup[,.N],       counts_up=cup[,count],       index_up=t(data.matrix(cup[,.(id1,id2)])), dist_up=cup[,distance],
               Ndown=cdown[,.N],   counts_down=cdown[,count],   index_down=t(data.matrix(cdown[,.(id1,id2)])), dist_down=cdown[,distance],
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag_centered=opt$par$beta_diag_centered)
  optimizing(model, data=data, as_vector=F, hessian=F, iter=1, verbose=verbose, init=0)
}

csnorm_predict_binned = function(model, biases, counts, opt, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL,
                               bf_per_decade=5, verbose=F, circularize=-1L) {
  csub=copy(counts) #need to implement taking only needed part of matrix, and reporting log_nu and log_delta appropriately
  bsub=copy(biases)
  Kdiag=round((log10(opt$dmax)-log10(opt$dmin))*bf_per_decade)
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
  data = list( Kdiag=Kdiag, npoints=npoints, circularize=circularize,
               S1=bsub[!is.na(bin1),.N], S2=bsub[!is.na(bin2),.N], 
               cutsites1=bsub[!is.na(bin1),pos], cutsites2=bsub[!is.na(bin2),pos],
               dmin=opt$dmin, dmax=opt$dmax,
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
    a=integrate(function(x){exp(dgamma(x,a2,rate=b2,log=T)+pgamma(x,a1,rate=b1,log.p=T))},xmin,xmax, stop.on.error=F)
    if (a$abs.error<=0 | a$message != "OK") {NA} else {a$value}
  }
}

#compute p(normal2>normal1) = \int_{-infty}^{+infty} dx p_normal2(x) \int_{-infty}^{x} dy p_normal1(y)
compute_normal_overlap = function(mu1,sd1,mu2,sd2, bounds=5, ncores=1) {
  registerDoParallel(cores=ncores)
  foreach (m1=mu1, s1=sd1, m2=mu2, s2=sd2, .packages="stats", .combine=c) %dopar% {
    #infinite integral does not work very well so truncate outer integral
    xmin = min(m1-bounds*s1,m2-bounds*s2)
    xmax = max(m1+bounds*s1,m2+bounds*s2)
    a=integrate(function(x){exp(dnorm(x,mean=(m2-m1)/s1,sd=s2/s1,log=T)+pnorm(x,mean=0,sd=1,log.p=T))},(xmin-m1)/s1,(xmax-m1)/s1)
    if (a$abs.error<=0 | a$message != "OK") {NA} else {a$value}
  }
}

detect_interactions = function(binned, dispersions, threshold=0.95, ncores=1, normal.approx=100){
  #report gamma parameters
  mat=copy(binned)
  mat[,c("alpha1","beta1"):=list(dispersions,dispersions/expected)]
  mat[,c("alpha2","beta2"):=list(alpha1+observed,beta1+1)]
  mat[alpha1<normal.approx,c("prob.observed.gt.expected","detection.type"):=list(compute_gamma_overlap(alpha1,beta1,alpha2,beta2,ncores=ncores),"gamma")]
  mat[,c("mean1","sd1"):=list(expected,expected/sqrt(dispersions))]
  mat[,c("mean2","sd2"):=list(alpha2/beta2, sqrt(alpha2)/beta2)]
  mat[alpha1>=normal.approx,c("prob.observed.gt.expected","detection.type"):=list(compute_normal_overlap(mean1,sd1, mean2, sd2, ncores=ncores),"normal")]
  mat[,prob.observed.gt.expected:=as.numeric(prob.observed.gt.expected)]
  mat[is.na(prob.observed.gt.expected)&observed>=expected,c("prob.observed.gt.expected","detection.type"):=list(ppois(observed,expected,lower.tail=F),"poisson")]
  mat[is.na(prob.observed.gt.expected)&observed<expected,c("prob.observed.gt.expected","detection.type"):=list(ppois(observed,expected,lower.tail=T),"poisson")]
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

read_and_prepare = function(infile, outprefix, skip=0L, both=T, distance_bins_per_decade=100, circularize=-1, dangling.L = c(0,4), dangling.R = c(3,-1)) {
  message("*** READ")
  data=read_tsv(infile, skip=skip)
  message("*** CATEGORIZE")
  data = categorize_by_new_type(data, dangling.L = dangling.L, dangling.R = dangling.R)
  message("*** BIASES AND COUNTS")
  cs_data = prepare_for_sparse_cs_norm(data, both=both, circularize=circularize)
  dset_statistics(cs_data$biases,cs_data$counts)
  message("*** WRITE")
  write.table(cs_data$biases, file = paste0(outprefix,"_biases.dat"), quote=F, row.names = F)
  write.table(cs_data$counts, file = paste0(outprefix,"_counts.dat"), quote=F, row.names = F)
  if (both==T) {
    both=cs_data$both
    save(both, file = paste0(outprefix,"_both.RData"))
  }
  return(cs_data)
}

postprocess = function(biases, counts, op, resolution=10000, ncores=30, predict.all.means=T, circularize=-1L) {
  smpred = stan_model(file = "cs_norm_predict.stan")
  smbin = stan_model("cs_norm_predict_binned.stan")
  smdisp = stan_model("cs_norm_binned_dispersions.stan")
  ### run remaining steps
  if (predict.all.means==T) {
    message("*** predict all means")
    op$pred=csnorm_predict_all(smpred, biases, counts, op, verbose=T)$par
  }
  message("*** buid binned matrices")
  op$binned=csnorm_predict_binned(smbin, biases, counts, op, resolution=resolution, circularize=circularize)
  message("*** estimate dispersions")
  op$disp=get_dispersions(smdisp, op$binned$mat)$par
  message("*** detect interactions")
  op$mat=detect_interactions(op$binned$mat, op$disp$dispersion, ncores=ncores) #interaction detection using binned dispersion estimates
  return(op)
}

diagnose_biases = function(outprefix, coverage=4, square.size=150000) {
  diagsquares = run_split_parallel_squares(biases, square.size, coverage, diag.only=T)$diagsquares
  p0=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2),alpha=0.1)+stat_function(fun=identity)
  diagsquares[,fname:=paste0(outprefix,"_biases_op_",begin1,".RData")]
  ops = sapply(X=diagsquares[,fname], FUN=function(x){ if (file.exists(x)) {a=load(x); return(get(a[1]));} else {return(NA)}}, USE.NAMES=T, simplify=F)
  diagsquares[,runtime:=sapply(X=ops,FUN=function(x){if (class(x)=="list") {return(x$runtime)} else {return(NA)}})]
  diagsquares[,result:=as.factor(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=1) {return(tail(x$output,1)[[1]])} else {return(NA)}} else {return(NA)}}))]
  diagsquares[,nsteps:=as.integer(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=3) {return(strsplit(tail(x$output,3)[[1]],'\\s+')[[1]][2])} else {return(NA)}} else {return(NA)}}))]
  p1=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=nsteps))+stat_function(fun=identity)
  p2=ggplot(diagsquares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=result))+stat_function(fun=identity)
  return(list(p0,p1,p2))
}

diagnose_counts = function(outprefix, coverage.extradiag=1, square.size=150000) {
  squares = run_split_parallel_squares(biases, square.size, coverage.extradiag, diag.only=F)$squares
  squares[,c("id1","id2"):=list(as.integer(factor(begin1)),as.integer(factor(begin2)))]
  p0=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2))+stat_function(fun=identity)
  #ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=is.na(nsteps)&begin2-begin1<1500000))+stat_function(fun=identity)
  squares[,fname:=paste0(outprefix,"_counts_op_",begin1,"_",begin2,".RData")]
  ops = sapply(X=squares[,fname], FUN=function(x){ if (file.exists(x)) {a=load(x); return(get(a[1]));} else {return(NA)}}, USE.NAMES=T, simplify=F)
  squares[,runtime:=sapply(X=ops,FUN=function(x){if (class(x)=="list") {return(x$runtime)} else {return(NA)}})]
  squares[,result:=as.factor(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=1) {return(tail(x$output,1)[[1]])} else {return(NA)}} else {return(NA)}}))]
  squares[,nsteps:=as.integer(sapply(X=ops,FUN=function(x){if (class(x)=="list") {
    if (length(x$output)>=3) {return(strsplit(tail(x$output,3)[[1]],'\\s+')[[1]][2])} else {return(NA)}} else {return(NA)}}))]
  p1=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=nsteps))+stat_function(fun=identity)
  p2=ggplot(squares)+geom_rect(aes(xmin=begin1,xmax=end1,ymin=begin2,ymax=end2,fill=result))+stat_function(fun=identity)
  return(list(p0,p1,p2))
}


#read_and_prepare("/scratch/ledily/mapped/HindIII_T0_both_filled_map_chr6.tsv", "data/ledily_T0_HindIII_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))
#read_and_prepare("/scratch/ledily/mapped/HindIII_T60_both_filled_map_chr6.tsv", "data/ledily_T60_HindIII_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))

#read_and_prepare("/scratch/ledily/mapped/NcoI_T0_both_filled_map_chr6.tsv", "data/ledily_T0_NcoI_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))
#read_and_prepare("/scratch/ledily/mapped/NcoI_T60_both_filled_map_chr6.tsv", "data/ledily_T60_NcoI_chr6", skip="SRR", dangling.L = c(0,4,5), dangling.R = c(3,-1,-2))

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

#prefix="ledily_T60_HindIII_chr6"
prefix="caulo_NcoI_all"
#prefix="ralph_Bcell_Sox2"
biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))

#ledily:  tad 821 151000000-153000000
biases=biases[pos>=151000000&pos<=153000000]
#biases=biases[pos>=150000000&pos<=154000000]

#ralph: restrict to 30500000,32500000
#biases=biases[pos<=30600714]
#biases=biases[pos>=30500000&pos<=32500000]
#biases=biases[pos>=30500000&pos<=30600000]
beginrange=biases[1,id]
endrange=biases[.N,id]
biases[,id:=id-beginrange+1]
prefix=biases[,paste0(prefix,"_",min(pos),"-",max(pos))]
counts=counts[id1>=beginrange&id1<=endrange&id2>=beginrange&id2<=endrange]
counts[,c("id1","id2"):=list(id1-beginrange+1,id2-beginrange+1)]

#caulo/toy: restrict to 100 consecutive rsites
biases=biases[pos>=1000000&pos<=3000000]
beginrange=biases[1,id]
endrange=biases[.N,id]
beginrange=1
endrange=biases[pos<4042929/2][,.N]
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


counts=fill_zeros(counts,biases)
counts[,distance:=abs(pos2-pos1)] #not filled
write.table(biases, file = paste0("data/",prefix,"_biases.dat"), quote=F, row.names = F)
write.table(counts, file = paste0("data/",prefix,"_counts.dat"), quote=F, row.names = F)

counts[,distance:=pmin(abs(pos2-pos1),4042929-abs(pos2-pos1)+1)] #not filled


### effect of parallelization
smfit = stan_model(file = "cs_norm_fit.stan")
opall <- csnorm_fit(smfit, biases, counts, bf_per_kb=1,
                                bf_per_decade=5, verbose = T, iter=100000)
opall=postprocess(biases, counts, opall, resolution=10000, ncores=30)
#save(opall, file = paste0("data/",prefix,"_op_maxcount_-1_redo.RData"))
load(paste0("data/",prefix,"_op_maxcount_-1.RData"), verbose=T)

coverage=4
square.size=150000
oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=1,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=100000, ncpus=30, homogenize=F)
oppar=postprocess(biases, counts, oppar, resolution=10000, ncores=30)
#save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb10.RData"))

load("data/",prefix,"_op_maxcount_-1.RData")
opserial=op
load(paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k.RData"))

#lFC histogram and matrices
ggplot(rbind(opserial$mat[,.(bin1,bin2,dset="serial",lFC)],oppar$mat[,.(bin1,bin2,dset="parallel",lFC)]))+geom_histogram(aes(lFC,fill=dset),position="dodge")
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_hist.png"), width=10, height=7.5)
ggplot()+geom_raster(data=oppar$mat,aes(begin1,begin2,fill=lFC))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+
  geom_raster(data=opserial$mat,aes(begin2,begin1,fill=lFC))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opserial$mat[is.interaction==T])
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_lFC_mat.png"), width=10, height=7.5)
#normalized matrix
ggplot()+geom_raster(data=oppar$mat,aes(begin1,begin2,fill=log(normalized)))+geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected>0.5),data=oppar$mat[is.interaction==T])+
  geom_raster(data=opserial$mat,aes(begin2,begin1,fill=log(normalized)))+geom_point(aes(begin2,begin1,colour=prob.observed.gt.expected>0.5),data=opserial$mat[is.interaction==T])
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_normalized.png"), width=10, height=7.5)
#genomic biases
#ggplot(melt(data.table(id=biases[,pos],serial=exp(opserial$par$eRJ+opserial$par$log_nu),parallel=exp(oppar$par$eRJ+oppar$par$log_nu)),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))+geom_point(data=biases,aes(pos,rejoined))
ggplot(melt(data.table(id=biases[,pos],serial=opserial$par$log_nu,parallel=oppar$par$log_nu),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu.png"), width=10, height=7.5)
ggplot(melt(data.table(id=biases[,pos],serial=opserial$par$log_delta,parallel=oppar$par$log_delta),id.vars="id",variable.name="delta"))+geom_line(aes(id,value,colour=delta))
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_delta.png"), width=10, height=7.5)
ggplot(rbind(data.table(bias="nu",serial=opserial$par$log_nu,parallel=oppar$par$log_nu),data.table(bias="delta",serial=opserial$par$log_delta,parallel=oppar$par$log_delta)))+
  geom_point(aes(serial,parallel))+facet_grid(.~bias)+stat_function(fun=identity)
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_nu_delta.png"), width=10, height=7.5)
#diagonal biases
ggplot(melt(data.table(dist=counts[,distance],serial=opserial$pred$log_decay_close,parallel=oppar$pred$log_decay_far, key="dist"),id.vars="dist"))+
  geom_line(aes(dist,value,colour=variable))+scale_x_log10()
ggsave(filename=paste0("images/",prefix,"_serial_vs_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_decay.png"), width=10, height=7.5)
c(opserial$mat[,mean(lFC)], oppar$mat[,mean(lFC)])
c(opserial$par$eRJ,oppar$par$eRJ)
c(opserial$par$eDE,oppar$par$eDE)
c(opserial$par$eC,oppar$par$eC)
c(opserial$disp$alpha,oppar$disp$alpha)




### single run plots
coverage=4
square.size=150000
bf_per_kb=1
circularize=4042929
oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=10000, ncpus=30,
                         homogenize=F, outprefix="tmp/test", circularize=circularize)
#oppar=run_split_parallel_recovery(counts, biases, outprefix, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
#                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=10000, ncpus=30, homogenize=F, circularize=circularize)
oppar=postprocess(biases, counts, oppar, resolution=10000, ncores=30, circularize=circularize, predict.all.means=F)
oppar$ice=iterative_normalization(oppar$mat, niterations=1, resolution=10000, return.binned=T)
save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))

load("data/caulo_NcoI_3189-2020271_op_maxcount_-1_parallel_inhomogeneous_cov4X_sq150k_bfpkb01.RData")
load("data/caulo_BglIIr1_1-2009521_op_maxcount_-1_parallel_inhomogeneous_cov4X_sq150k_bfpkb01.RData")
load("data/caulo_BglIIr2_15389-2009521_op_maxcount_-1_parallel_inhomogeneous_cov4X_sq150k_bfpkb01.RData")

#prefix="ledily_T0_NcoI_chr6_150000151-153992077"
#prefix="ledily_T0_NcoI_chr6_151001273-152989574"
#prefix="ledily_T0_HindIII_chr6_150004969-153999619"
#prefix="ledily_T0_HindIII_chr6_151000912-152997061"
#prefix="ledily_T60_NcoI_chr6_150000151-153992077"
#prefix="ledily_T60_NcoI_chr6_151001273-152989574"
#prefix="ledily_T60_HindIII_chr6_150004969-153999619"
#prefix="ledily_T60_HindIII_chr6_151000912-152997061"
#prefix="caulo_BglIIr1_1-2009521"
#prefix="caulo_BglIIr2_15389-2009521"
#prefix="caulo_NcoI_3189-2020271"
prefix="caulo_NcoI_all"
#prefix="caulo_NcoI_1004680-2998140"

biases=fread(paste0("data/",prefix,"_biases.dat"))
counts=fread(paste0("data/",prefix,"_counts.dat"))
load(paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))
setkey(oppar$mat, begin1,begin2)
#matrix and detected interactions
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(normalized)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(normalized)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none")
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_normalized.png"), width=10, height=9)
#lFC
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=lFC))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=lFC))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_lFC_mat.png"), width=10, height=9)
#raw
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(observed)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+scale_fill_gradient(low="white", high="black")+theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_observed.png"), width=10, height=9)
#ice
ggplot()+geom_raster(data=oppar$ice, aes(begin1,begin2,fill=log(N)))+geom_raster(data=oppar$ice, aes(begin2,begin1,fill=log(N)))+scale_fill_gradient(low="white", high="black")+  theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_ICE.png"), width=10, height=9)+  theme(legend.position = "none") 
#ice vs raw
ggplot()+geom_raster(data=oppar$ice, aes(begin1,begin2,fill=log(N)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+scale_fill_gradient(low="white", high="black")+  theme(legend.position = "none") 
ggsave(filename = paste0("images/",prefix,"_ICE_vs_observed.png"), width=10, height=9)+  theme(legend.position = "none") 
#cs vs ice
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(normalized)))+geom_raster(data=oppar$ice, aes(begin2,begin1,fill=log(N)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none") 
#observed vs expected
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(expected)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(observed)))+
  #geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+
  scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none") 
#nu
ggplot(melt(data.table(id=biases[,pos],parallel=oppar$par$log_nu),id.vars="id",variable.name="nu"))+geom_line(aes(id,value,colour=nu))
#decay
ggplot(data.table(dist=oppar$binned$distance,decay=oppar$binned$decay)[dist>1e3])+geom_line(aes(dist,decay))+scale_x_log10()+scale_y_log10()
#decay and normalized counts
#oppar$mat[,distance:=pmin(begin2-begin1,4042929+1-(begin2-begin1))]
#setkey(oppar$mat,distance)
#ggplot(data.table(dist=oppar$binned$distance,decay=oppar$binned$decay)[dist>1e3])+geom_line(aes(dist,decay),colour="red")+scale_x_log10()+scale_y_log10()+
#  geom_point(data=oppar$mat,aes(distance,normalized),alpha=0.01)
ggplot(data.table(normalized=counts[,contact.close]*exp(oppar$pred$log_decay_close-oppar$pred$log_mean_cclose),distance=counts[,distance],
                  decay=exp(oppar$pred$log_decay_close))[distance>1000][sample(.N,min(.N,100000))])+
  geom_point(aes(distance,normalized),alpha=0.01)+
  geom_line(aes(distance,decay),colour="red")+scale_x_log10()+scale_y_log10()
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_decay_with_counts.png"), width=10, height=9)
#lfC hist
ggplot(oppar$mat)+geom_histogram(aes(lFC))
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_lFC_hist.png"), width=10, height=9)
#detected interactions and 95% CI for posterior
oppar$mat[,upper:=qgamma(0.975,shape=alpha1,rate=beta1)]
oppar$mat[,lower:=qgamma(0.025,shape=alpha1,rate=beta1)]
oppar$mat[,CI:=upper-lower]
ggplot()+geom_raster(data=oppar$mat, aes(begin1,begin2,fill=log(normalized)))+geom_raster(data=oppar$mat, aes(begin2,begin1,fill=log(CI*normalized/observed)))+
  geom_point(aes(begin1,begin2,colour=prob.observed.gt.expected<0.5),data=oppar$mat[is.interaction==T])+scale_fill_gradient(low="white", high="black")+
  theme(legend.position = "none")#+scale_colour_brewer(type="div",palette="RdBu")
ggsave(filename = paste0("images/",prefix,"_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_normalized_with_error.png"), width=10, height=9)



### reload intermediate files of long run
#plot bias runs
ops=sapply(X=Sys.glob("ralph_*biases_op_*.RData"), FUN=function(x){ a=load(x); return(get(a[1]))}, USE.NAMES=T, simplify=F)
failures=data.table(failed=sapply(ops, function(x){length(grep("failed", tail(x$output,1)))>0}),
                    runtime=sapply(ops, function(x){x$runtime}),
                    sqsizes=tstrsplit(names(ops),"_")[[7]],
                    dispersion=sapply(ops, function(x){x$par$alpha}),
                    deviance=sapply(ops, function(x){x$par$deviance_proportion_explained}))
failures[,sqsizes:=as.numeric(sapply(sqsizes, function(x){substr(x,1,nchar(x)-6)}))]
ggplot(failures)+geom_point(aes(sqsizes, runtime, colour=failed))
ggplot(failures)+geom_point(aes(sqsizes, deviance, colour=failed))
ggplot(failures)+geom_point(aes(sqsizes, dispersion, colour=failed))+scale_y_log10()


#plot counts runs
ops=sapply(X=Sys.glob("tmp/ralph_*counts_op_*.RData"), FUN=function(x){ a=load(x); return(get(a[1]))}, USE.NAMES=T, simplify=F)
failures=data.table(failed=sapply(ops, function(x){length(grep("failed", tail(x$output,1)))>0}),
                    runtime=sapply(ops, function(x){x$runtime}),
                    sqsizes=tstrsplit(names(ops),"_")[[7]],
                    deviance=sapply(ops, function(x){x$par$deviance_proportion_explained}))
failures[,sqsizes:=as.numeric(sapply(sqsizes, function(x){substr(x,1,nchar(x)-6)}))]
failures[,idx:=.I]
ggplot(failures)+geom_point(aes(idx, runtime, colour=failed))
ggplot(failures)+geom_point(aes(idx, deviance, colour=failed))

ops.bias=sapply(X=Sys.glob(paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_biases_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
ops.bias = output.binder(ops.bias)
save(ops.bias, file=paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_ops_bias.RData"))

ops.count=sapply(X=Sys.glob(paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_counts_ret_*.RData")), FUN=function(x){ a=load(x); return(list(ret=get(a[1]), out="blah", runtime=-1))}, USE.NAMES=T, simplify=F)
ops.count = output.binder(ops.count)
save(ops.count, file=paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,"_ops_count.RData"))

oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=100000, ncpus=30, homogenize=F, ops.count=ops.count, ops.bias=ops.bias)
oppar=postprocess(biases, counts, oppar, resolution=50000, ncores=30, predict.all.means=T)
oppar$ice=iterative_normalization(oppar$mat, niterations=1, resolution=50000, return.binned=T)
save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))


### generate LGF matrix
a=fread("~/simulations/caulobacter/data/GSM1120448_Laublab_NcoI_HiC_NA1000_swarmer_cell_untreated_overlap_after_HiCNorm_fullgen.count")
setnames(a,c("bin1","bin2","lgfcount"))
bins=data.table(begin=seq(1,4045000+10000,10000))
bins[,name:=paste0("bin",.I-1)]
setkey(bins,name)
setkey(a,bin1)
a=bins[a]
setnames(a,c("begin", "name"), c("begin1","bin1"))
setkey(a,bin2)
a=bins[a]
setnames(a,c("begin", "name"), c("begin2","bin2"))
#newbins=seq(1,4045000+20000,20000)
#a[,nb1:=cut2(begin1,newbins)]
#a[,nb2:=cut2(begin2,newbins)]
#a=a[,.(min(begin1),min(begin2),lgfcount=sum(lgfcount)),by=c("nb1","nb2")]
#a[,c("ign","begin1","ign2"):=tstrsplit(nb1,'[[,]')]
#a[,c("ign","begin2","ign2"):=tstrsplit(nb2,'[[,]')]
#a[,c("begin1","begin2"):=list(as.numeric(begin1),as.numeric(begin2))]
#ggplot(a[begin1<=2003200&begin2<=2003200])+geom_raster(aes(begin1,begin2,fill=log(lgfcount)))+geom_raster(aes(begin2,begin1,fill=log(lgfcount)))+
#  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")
ggplot(a)+geom_raster(aes(begin1,begin2,fill=log(lgfcount)))+geom_raster(aes(begin2,begin1,fill=log(lgfcount)))+
  scale_fill_gradient(low="white", high="black")+theme(legend.position = "none")
ggsave(filename = paste0("images/",prefix,"_HiCNorm.png"), width=10, height=9)+  theme(legend.position = "none") 
