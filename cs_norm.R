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
registerDoParallel(cores=30)

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

iterative_normalization = function(data, niterations=100, resolution=1000000, return.binned=F) {
  #bin data and remove diagonal
  genome.len=64444167
  bins=seq(0,genome.len+resolution,resolution)
  bdata = data[category %in% c("contact close","contact up","contact down","contact far"),
               .(begin1,begin2,bin1=cut2(begin1, bins, oneval=F, onlycuts=T, digits=10),
                 bin2=cut2(begin2, bins, oneval=F, onlycuts=T, digits=10))]
  binned = bdata[,.N,by=c("bin1","bin2")]
  binned = binned[bin1!=bin2]
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

stan_matrix_to_datatable = function(opt, x) {
  vals=data.table(opt)
  vals[,x:=x]
  melt(data.table(vals), id.vars="x")
}
optimize_all_meanfield = function(model, biases, counts, meanfield, maxcount, bf_per_kb=1, bf_per_decade=5,
                                  iter=10000, verbose=T, mincount=-1, ...) {
  cclose=counts[contact.close>maxcount,.(id1,id2,distance,count=contact.close)][count>mincount]
  cfar=counts[contact.far>maxcount,.(id1,id2,distance,count=contact.far)][count>mincount]
  cup=counts[contact.up>maxcount,.(id1,id2,distance,count=contact.up)][count>mincount]
  cdown=counts[contact.down>maxcount,.(id1,id2,distance,count=contact.down)][count>mincount]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round(counts[,(log10(max(pos2-pos1))-log10(min(pos2-pos1)))*bf_per_decade])
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag,
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
  optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, ...)
}


vb_all_meanfield = function(model, biases, counts, meanfield, maxcount, bf_per_kb=1, bf_per_decade=5,
                                  iter=10000, mincount=-1, ...) {
  cclose=counts[contact.close>maxcount,.(id1,id2,distance,count=contact.close)][count>mincount]
  cfar=counts[contact.far>maxcount,.(id1,id2,distance,count=contact.far)][count>mincount]
  cup=counts[contact.up>maxcount,.(id1,id2,distance,count=contact.up)][count>mincount]
  cdown=counts[contact.down>maxcount,.(id1,id2,distance,count=contact.down)][count>mincount]
  mf=list()
  mf$Nkl=meanfield$Nkl[count<=maxcount]
  mf$Nkr=meanfield$Nkr[count<=maxcount]
  mf$Nkd=meanfield$Nkd[count<=maxcount]
  Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
  Kdiag=round(counts[,(log10(max(pos2-pos1))-log10(min(pos2-pos1)))*bf_per_decade])
  data = list( Krow=Krow, S=biases[,.N],
               cutsites=biases[,pos], rejoined=biases[,rejoined],
               danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
               Kdiag=Kdiag,
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
  op=vb(model, data=data, iter=iter, ...)
  list(par=sapply(op@model_pars, function(x){as.numeric(get_posterior_mean(op, x))}))
}





predict_full = function(model, biases, counts, opt, verbose=T) {
  #need to ensure that counts span exactly the same distance range as during fitting
  data = list( Kdiag=length(op$par$beta_diag)+1, S=biases[,.N], cutsites=biases[,pos], N=counts[,.N],
               counts=t(data.matrix(counts[,.(contact.close,contact.far,contact.up,contact.down)])),
               cidx=t(data.matrix(counts[,.(id1,id2)])),
               eC=opt$par$eC, log_nu=opt$par$log_nu, log_delta=opt$par$log_delta,
               beta_diag=opt$par$beta_diag, alpha=opt$par$alpha)
  optimizing(model, data = data, as_vector=F, hessian=F, iter=1, verbose=verbose)
}


bin_counts = function(counts, biases, resolution, b1=NULL, b2=NULL, e1=NULL, e2=NULL, normalized=T) {
  if (is.null(b1)) b1=counts[,min(pos1)]-1
  if (is.null(b2)) b2=counts[,min(pos2)]-1
  if (is.null(e1)) e1=counts[,max(pos1)]+1
  if (is.null(e2)) e2=counts[,max(pos2)]+1
  bins1=seq(b1,e1+resolution,resolution)
  bins2=seq(b2,e2+resolution,resolution)
  #
  counts.wrapped=rbindlist(list(counts[,.(pos1,pos2,contact.close,log_decay,log_mean_cclose)],
                                counts[,.(pos1,pos2,contact.far,log_decay,log_mean_cfar)],
                                counts[,.(pos1,pos2,contact.up,log_decay,log_mean_cup)],
                                counts[,.(pos1,pos2,contact.down,log_decay,log_mean_cdown)]))
  setnames(counts.wrapped, c("pos1","pos2","count","log_decay","log_mean"))
  if (normalized==T) {
    sub = counts.wrapped[,.(pos1,pos2,bin1=cut2(pos1, bins1, oneval=F, onlycuts=T, digits=10, minmax=F),
                            bin2=cut2(pos2,bins2, oneval=F, onlycuts=T, digits=10, minmax=F),
                            weight=exp(log_mean-log_decay))
                         ][,.(N=sum(weight, na.rm = T)),by=c("bin1","bin2")]
  } else {
    sub = counts.wrapped[,.(pos1,pos2,bin1=cut2(pos1, bins1, oneval=F, onlycuts=T, digits=10, minmax=F),
                            bin2=cut2(pos2,bins2, oneval=F, onlycuts=T, digits=10, minmax=F),
                            weight=count)
                         ][,.(N=sum(weight, na.rm = T)),by=c("bin1","bin2")]
  }
  #remove NAs in case a zoom was performed
  sub=sub[complete.cases(sub)]
  #divide by number of rsites in each bin
  if (normalized==T) {
    ns1 = biases[,.(bin1=cut2(pos, bins1, oneval=F, onlycuts=F, digits=10))][,.(nsites1=.N),keyby=bin1]
    setkey(sub, bin1)
    sub=ns1[sub]
    ns2 = biases[,.(bin2=cut2(pos, bins2, oneval=F, onlycuts=F, digits=10))][,.(nsites2=.N),keyby=bin2]
    setkey(sub, bin2)
    sub=ns2[sub]
    sub[,N:=N/(nsites1*nsites2)]
  }
  #write begins/ends
  bin1.begin=sub[,bin1]
  bin1.end=sub[,bin1]
  bin2.begin=sub[,bin2]
  bin2.end=sub[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[[,]")[2][[1]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[[,)]")[2][[1]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[[,]")[2][[1]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[[,)]")[2][[1]]
  sub[,begin1:=as.integer(as.character(bin1.begin))]
  sub[,end1:=as.integer(as.character(bin1.end))]
  sub[,begin2:=as.integer(as.character(bin2.begin))]
  sub[,end2:=as.integer(as.character(bin2.end))]
  return(sub)
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
  mcounts[,bdist:=cut(distance,dbins,ordered_result=T,right=F,include.lowest=T)]
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

todo = function(){
  gammaQuery=function(mu,theta,count){
    alpha1=theta
    alpha2=theta+count
    beta1=theta/mu
    beta2=beta1+1
    return(function(x){dgamma(x,alpha2,rate=beta2)*pgamma(x,alpha1,rate=beta1)})
  }
  igammaQuery_c=function(mu,theta,threshold=0.95){
    function(x){integrate(gammaQuery(mu=mu,theta=theta,count=x),0,+Inf)$value-threshold}
  }
  #uniroot(igammaQuery_c(mu=2,theta=2,threshold=0.9),c(0,10))$root
  #uniroot(igammaQuery_c(mu=2,theta=20,threshold=0.9),c(10,30))$root
  igammaQuery_theta=function(mu,count,threshold=0.95){
    function(x){integrate(gammaQuery(mu=mu,theta=x,count=count),0,+Inf)$value-threshold}
  }
  #uniroot(igammaQuery_theta(mu=2,count=3,threshold=0.9),c(0,30))$root
  #uniroot(igammaQuery_theta(mu=2,count=10,threshold=0.9),c(1,30))$root
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

#read_and_prepare("/scratch/rao/mapped/HICall_both_filled_map_chr1.tsv", "data/rao_HICall_chr1_all", skip="SRR")

read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate1_reads_int.tsv",
                 "data/caulo_BglIIr1_all", skip="SRR", circularize=4042929)
read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_BglII_replicate2_reads_int.tsv",
                 "data/caulo_BglIIr2_all", skip="SRR", circularize=4042929)
read_and_prepare("/scratch/caulobacter/6_preprocessing_raw_reads/3_InteractionMaps/Caulobacter_NcoI_reads_int.tsv",
                 "data/caulo_NcoI_all", skip="SRR", circularize=4042929)

dset=generate_fake_dataset(num_rsites=100, genome_size = 100000, eC = .1, eRJ = 5, eDE = 7)
counts=dset$counts
biases=dset$biases
dset_statistics(biases,counts)
meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)


prefix="caulo_NcoI_all"
biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))
#meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 100)
#save(meanfield, file = paste0("data/",prefix,"_meanfield_100.RData"))
load(paste0("data/",prefix,"_meanfield_100.RData"))

#restrict to 100 consecutive rsites
beginrange=150
endrange=170
biases=biases[beginrange:endrange]
biases[,id:=id-beginrange+1]
prefix=biases[,paste0("caulo_NcoI_",min(pos),"-",max(pos))]
counts=counts[id1>=beginrange&id1<=endrange&id2>=beginrange&id2<=endrange]
counts[,c("id1","id2"):=list(id1-beginrange+1,id2-beginrange+1)]
meanfield$Nkl=meanfield$Nkl[id>=beginrange&id<=endrange][,.(id=id-beginrange+1,count,N)]
meanfield$Nkr=meanfield$Nkr[id>=beginrange&id<=endrange][,.(id=id-beginrange+1,count,N)]

#### optimization wihout prior guesses
smfit = stan_model(file = "cs_norm_fit.stan")
smpred = stan_model(file = "cs_norm_predict.stan")
maxcount=4
a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                         bf_per_decade=5, verbose = T, iter=100000, tol_rel_grad=1e3, tol_rel_obj=1e3))
op.slice <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                           bf_per_decade=5, verbose = T, iter=100000, tol_rel_grad=1e3, tol_rel_obj=1e3)
a=data.table(logmean.redo=op$par$log_mean_cclose,
             logmean.slice=op.slice$par$log_mean_cclose)
ggplot(a)+geom_point(aes(logmean.redo,logmean.slice))

a=melt(data.table(lognu.short1=op.short1$par$eC+op.short1$par$log_nu,
                  lognu.short2=op.short1$par$eC+op.short2$par$log_nu,
                  lognu.long1=op.long1$par$eC+op.long1$par$log_nu,
                  lognu.long2=op.long2$par$eC+op.long2$par$log_nu,
                  bias=150:170), id.vars="bias")
ggplot(a)+geom_line(aes(bias,value,colour=variable))


b=system.time(op <- vb_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                           bf_per_decade=5, iter=100000, tol_rel_obj=1e3, algorithm="fullrank"))
op.slice <- vb_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                     bf_per_decade=5, iter=100000, tol_rel_obj=1e3, algorithm="fullrank")


#op$par$time=a[1]+a[4]
#op$par$maxcount=maxcount
#op.pred <- predict_full(smpred, biases, counts, op, verbose = T)

a=melt(data.table(lognu.redo=op$par$eC+op$par$log_nu,
                  lognu.slice=op.slice$par$eC+op.slice$par$log_nu,
                  lognu.full=op.full$par$eC+op.full$par$log_nu[140:160],bias=140:160), id.vars="bias")
ggplot(a)+geom_line(aes(bias,value,colour=variable))

b=melt(data.table(lognu.redo=op$par$log_nu,
                  lognu.slice=op.slice$par$log_nu,
                  bias=140:160), id.vars="bias")
b[,sum(value),by=variable]
ggplot(b)+geom_line(aes(bias,value,colour=variable))



#save(op, file=paste0("data/",prefix,"_op_maxcount_",maxcount,".RData"))

meanfield=bin_for_mean_field(biases, counts, distance_bins_per_decade = 10)
smfit = stan_model(file = "sparse_cs_norm_fit_meanfield.stan")
smpred = stan_model(file = "sparse_cs_norm_predict.stan")
ops = foreach (maxcount=c(-1,0,1,2,3,4,6,1e18)) %:% foreach (repetition=1:10) %dopar% {
  message("***** ",maxcount)
  #fit
  a=system.time(op <- optimize_all_meanfield(smfit, biases, counts, meanfield, maxcount=maxcount, bf_per_kb=1,
                                             bf_per_decade=5, verbose = T, iter=10))
  op$par$time=a[1]+a[4]
  if (maxcount == -1) {
    op$par$maxcount="exact"
  } else if (maxcount>counts[,max(c(contact.close,contact.far,contact.down,contact.up))]) {
    op$par$maxcount="meanfield"
  } else {
    op$par$maxcount=maxcount
  }
  op$par$repetition=repetition
  #predict
  op.pred <- predict_full(smpred, biases, counts, op, verbose = T)
  op$pred=op.pred$par
  op#ops[[paste0("op_",maxcount,"_",repetition)]]=op
}
ops = unlist(ops, recursive=F)
ops=ops[!sapply(ops, is.null)] 
names(ops) <- seq(length(ops))
save(ops, file = "data/caulo_mf_3-3.2M.RData")
#load("data/caulo_mf_3-3.2M.RData")










#visualize binned matrices
binned = bin_counts(counts=counts, biases=biases, resolution=5000, normalized=F)
ggplot(binned, aes(begin1,begin2, fill=log(N)))+geom_raster()+scale_fill_gradient(low="white", high="black")
ggsave(filename = "images/rao_HICall_chrX_73780165-74230165_5k_raw.png", width=10, height=7.5)
write.table(binned[,.(begin1,end1,begin2,end2,N)], file = "data/binned_Rao_MboI_chrX_73780165-74230165_raw_5k.dat", quote = F, row.names = F)
#
binned = bin_counts(counts=counts, biases=biases, resolution=5000, normalized=T)
ggplot(binned, aes(begin1,begin2, fill=log(N)))+geom_raster()+scale_fill_gradient(low="white", high="black")
ggsave(filename = "images/rao_HICall_chrX_73780165-74230165_5k_div_nsites.png", width=10, height=7.5)
write.table(binned[,.(begin1,end1,begin2,end2,N)], file = "data/binned_Rao_MboI_chrX_73780165-74230165_div_nsites.dat", quote = F, row.names = F)




#precision and accuracy on single-valued params
single_params = data.table(
  sapply(c("eC","eRJ","eDE","lambda_nu","lambda_delta","lambda_diag","alpha","deviance_proportion_explained",
           "repetition","time"),
         function(y){sapply(names(ops), function(x){ops[[x]]$par[[y]]})}))
single_params[,maxcount:=sapply(names(ops), function(x){ops[[x]]$par$maxcount})]
single_params[,maxcount:=ordered(maxcount,levels=unique(maxcount))]
#single_params=single_params[eC<30&eRJ<6&eDE<7.5&lambda_delta<200] #caulo
#[maxcount!="meanfield"&eDE<20]
ggplot(melt(single_params[maxcount!="meanfield"],id.vars=c("maxcount","repetition")))+
  geom_boxplot(aes(maxcount,value,colour=(maxcount=="exact")))+
  facet_wrap(~variable, scales="free_y",shrink=T)+guides(colour=F)+
  #labs(title="Toy dataset, 40% zeros", x="mean field threshold", y=NULL)
  #labs(title="Caulobacter: 200k stretch, 40% zeros", x="mean field threshold", y=NULL)
  labs(title="Rao HICall chr20 300k stretch, 90% zeros", x="mean field threshold", y=NULL)
ggsave(filename = "images/rao_HICall_chr20_300k_mf_singleparams.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-3.2M_mf_singleparams.png", width=10, height=7.5)
ggsave(filename = "images/toy_high_mf_singleparams_lambda_decoupled_alpha_double_dbins_10_bfkb_1.png", width=10, height=7.5)

#precision on vector-valued params
all_params = data.table(
  sapply(names(ops[["1"]]$par), function(y){lapply(names(ops), function(x){ops[[x]]$par[[y]]})}))
all_params[,maxcount:=ordered(maxcount,levels=unique(maxcount))]
all_params[,repetition:=as.integer(repetition)]
mult_params=melt(all_params[,.(maxcount,repetition,log_nu,log_delta,beta_diag)], id.vars=c("maxcount","repetition"))
mult_params=mult_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("maxcount","repetition","variable")][
  ,.(std=sd(value)),by=c("maxcount","idx","variable")]
ggplot(mult_params)+
  geom_boxplot(aes(maxcount,std,colour=(maxcount=="exact")))+facet_wrap(~variable)+
  guides(colour=F)+scale_y_log10()+
  #labs(title="Toy dataset, 40% zeros",
  #labs(title="Caulobacter 200k stretch, 40% zeros: precision",
  labs(title="Rao HICall chr20 300k stretch, 90% zeros: precision",
       x="mean field threshold", y="standard deviation across repeats, for each coefficient")
ggsave(filename = "images/rao_HICall_chr20_300k_mf_multiparams_precision.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-3.2M_mf_multiparams_precision.png", width=10, height=7.5)
ggsave(filename = "images/toy_high_mf_multiparams_precision.png", width=10, height=7.5)

#accuracy on vector-valued params
ref_params=melt(all_params[maxcount=="exact",.(repetition,log_nu,log_delta,beta_diag)], id.vars="repetition")
ref_params=ref_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","repetition")][
  ,.(med=median(value)),keyby=c("idx","variable")]
mult_params=melt(all_params[,.(maxcount,repetition,log_nu,log_delta,beta_diag)], id.vars=c("maxcount","repetition"))
mult_params=mult_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("maxcount","repetition","variable")]
setkey(mult_params, idx, variable)
mult_params=ref_params[mult_params][,.(prec=dist(rbind(value,med))),by=c("maxcount","idx","variable")]
ggplot(mult_params)+
  geom_boxplot(aes(maxcount,prec,colour=(maxcount=="exact")))+facet_wrap(~variable)+guides(colour=F)+
  #labs(title="Toy dataset, 40% zeros",
  #labs(title="Caulobacter 200k stretch, 40% zeros: accuracy",
  labs(title="Rao HICall chr20 300k stretch, 90% zeros: accuracy",
       x="mean field threshold", y="distance from median of exact calculation")+coord_cartesian(ylim=c(0,10))
ggsave(filename = "images/rao_HICall_chr20_300k_mf_multiparams_accuracy.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-3.2M_mf_multiparams_accuracy.png", width=10, height=7.5)
ggsave(filename = "images/toy_high_mf_multiparams_accuracy.png", width=10, height=7.5)

#plot nu
pbegin=35100000 
pend=35200000 
#pbegin=3166716
#pend=3266716
#pbegin=3110583
#pend=3210583
#pbegin=0
#pend=20000
mean_params=melt(all_params[,.(maxcount,repetition,log_nu,log_delta,eRJ,eDE)], id.vars=c("maxcount","repetition"))
mean_params=mean_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","maxcount","repetition")][
  ,.(med=median(value)),keyby=c("maxcount","variable","idx")]
setkey(mean_params, idx)
ggplot(biases[pos>=pbegin&pos<=pend])+scale_y_log10()+coord_cartesian(ylim=c(0.1,10))+
  geom_point(aes(pos, dangling.L/exp(mean_params[variable=='eDE'&maxcount=="exact",med])),colour="orange")+
  geom_point(aes(pos, dangling.R/exp(mean_params[variable=='eDE'&maxcount=="exact",med])),colour="pink")+
  geom_point(aes(pos, rejoined/exp(mean_params[variable=='eRJ'&maxcount=="exact",med])),colour="red")+
  geom_line(data=biases[mean_params[variable=='log_nu']][pos>=pbegin&pos<=pend],aes(pos, exp(med), colour=maxcount))+
  scale_colour_brewer(palette="YlOrRd")+
  #geom_line(aes(pos, exp(true_log_nu)), colour="black")+
  #labs(title="Toy dataset, 40% zeros", x="genomic position", y="biases")
  labs(title="Rao HICall chr20 300k stretch: 100k example for nu", x="genomic position", y="biases")
#labs(title="Caulobacter: 100k example for nu", x="genomic position", y="biases")
ggsave(filename = "images/rao_HICall_chr20_300k_mf_nu_example.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-3.2M_mf_nu_example.png", width=10, height=7.5)
ggsave(filename = "images/toy_high_mf_nu_example.png", width=10, height=7.5)

#plot delta
ggplot(biases[pos>=pbegin&pos<=pend])+scale_y_log10()+
  geom_point(aes(pos, dangling.L/exp(mean_params[variable=='eDE'&maxcount=="exact",med])),colour="orange")+
  geom_line(data=dcast(biases[mean_params], ...~variable, value.var="med")[pos>=pbegin&pos<=pend],
            aes(pos, exp(log_nu+log_delta), colour=maxcount))+
  scale_colour_brewer(palette="YlOrRd")+
  #labs(title="Toy dataset, 40% zeros", x="genomic position", y="biases")
  labs(title="Rao HICall chr20 300k stretch: 100k example for delta", x="genomic position", y="biases")
#labs(title="Caulobacter: 100k example for delta", x="genomic position", y="biases")
ggsave(filename = "images/rao_HICall_chr20_300k_mf_delta_example.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-3.2M_mf_delta_example.png", width=10, height=7.5)
ggsave(filename = "images/toy_high_mf_delta_example.png", width=10, height=7.5)

#plot decay
decay_params = cbind(
  data.table(sapply(c("maxcount","repetition"), function(y){sapply(names(ops), function(x){ops[[x]]$par[[y]]})})),
  data.table(sapply(names(ops[["1"]]$pred), function(y){lapply(names(ops), function(x){ops[[x]]$pred[[y]][1:counts[id1==1,.N]]})})))
decay_params[,maxcount:=as.character(maxcount)]
decay_params[maxcount=="-1",maxcount:="exact"]
decay_params[,maxcount:=ordered(maxcount,levels=unique(maxcount))]
decay_params[,repetition:=as.integer(repetition)]
decay_params=melt(decay_params, id.vars=c("maxcount","repetition"))
decay_params=decay_params[,.(idx=c(1:length(unlist(value))),value=unlist(value)),by=c("variable","repetition","maxcount")]
decay_params=decay_params[,.(med=median(value)),keyby=c("idx","variable","maxcount")]
setkey(decay_params, idx)
decay_params[,dist:=rep(counts[id1==1,pos2-pos1],each=decay_params[idx==1,.N])]
decay_params=dcast(decay_params, ...~variable, value.var="med")
ggplot(decay_params)+scale_y_log10()+scale_x_log10()+ #[maxcount%in%c("exact","0","1")]
  scale_colour_brewer(palette="YlOrRd")+
  geom_line(aes(abs(dist), exp(log_decay), colour=maxcount))+
  #geom_line(data=counts[,.(decay=exp(true_log_decay)),by=distance],aes(distance,decay),colour="black")+
  labs(title="Rao HICall chr20 300k stretch: diagonal decay", x="distance", y="decay")
#labs(title="Caulobacter: 100k example for decay", x="genomic position", y="biases")
#labs(title="Toy dataset, high coverage, 40% zeros: decay", x="genomic position", y="biases")
ggsave(filename = "images/rao_HICall_chr20_300k_mf_decay.png", width=10, height=7.5)
ggsave(filename = "images/caulo_3-3.2M_mf_decay.png", width=10, height=7.5)
#

