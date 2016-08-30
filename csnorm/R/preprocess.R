#' @include csnorm.R
NULL

#' Read TADBit tsv file and return (paired-end) reads as data.table
#'
#' @param fname The filename
#' @param nrows,skip see \code{\link[data.table]{fread}}
#'
#' @return a data.table
#' @export
#'
#' @examples
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

#' Get subset of reads
#'
#' @param data 
#' @param b1,e1 begin / end of read range you want
#' @param b2,e2 (optional) if provided, give extradiagonal portion 
#'
#' @return data.table
#' @keywords internal
#' @export
#'
#' @examples
get_subset = function(data, b1, e1, b2=NULL, e2=NULL) {
  if (is.null(b2)) b2=b1
  if (is.null(e2)) e2=e1
  stopifnot(e1>b1)
  stopifnot(e2>b2)
  stopifnot(b2>=b1)
  return(data[begin1>=b1 & begin1+length1<=e1 & begin2>=b2 & begin2+length2<=e2])
}

#' Add category label to reads
#' 
#' Performs classification according to cut-site normalization scheme
#' 
#' @param sub reads data.table to cateogorize
#' @param maxlen positive integer.  Maximum admissible size of sonication 
#'   fragment. Fragments above this threshold will be categorized as "other".
#' @param dangling.L,dangling.R vectors of integers. Offset from cut site for
#'   reads to be considered dangling ends (Left or Right respectively)
#'   
#' @return The same data.table, with an additional "category"column
#' @export
#' 
#' @examples
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

#' Parse reads data.table and extract necessary objects for cut-site
#' normalization
#' 
#' @param data reads data.table
#' @param both boolean. Whether a merged table providing biases and their counts
#'   should be produced (default is FALSE as it is not used by cs norm)
#' @param circularize integer. Set this to the size of the chromosome if it is
#'   circular, otherwise leave as-is (default is -1)
#'   
#' @return A list containing biases and counts, two data tables representing the
#'   total number of reads of different categories at each cut site.
#' @export
#' 
#' @examples
prepare_for_sparse_cs_norm = function(data, both=F, circularize=-1) {
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
  X=rbind(X,rsites[!(id%in%X[,id]),.(id,pos=re.pos,dangling.L=0,dangling.R=0,rejoined=0)])
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

#' Print general statistics on the output of \code{\link{prepare_for_sparse_cs_norm}}
#'
#' @param biases,counts 
#'
#' @export
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

#' Generate a fake dataset that follows the cut site normalization model's specifications
#'
#' @param num_rsites Number of restriction sites
#' @param genome_size Size of the genome
#' @param eC,eRJ,eDE Exposure for counts, rejoined and dangling ends
#' @param alpha Dispersion
#'
#' @return a list of biases and counts, similar to \code{\link{prepare_for_sparse_cs_norm}}
#' @export
#'
#' @examples
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

#' Diagnostic plots to determine dangling.L, dangling.R and maxlen arguments
#'
#' @param infile A tadbit .tsv file
#' @param skip see \code{\link[data.table]{fread}}
#' @param nrows number of rows to read, default all.
#' @param window how many bases to plot around cut site
#' @param maxlen how many bases to plot off-diagonal
#'
#' @return Three plots, pleft and pright to determine dangling.L and dangling.R, and pdiag for maxlen
#' @export
#'
#' @examples
examine_dataset = function(infile, skip=0L, nrows=-1L, window=15, maxlen=1000) {
  data=read_tsv(infile, skip=skip, nrows=nrows)
  pleft=ggplot(data[abs(rbegin1-re.closest1)<=window])+geom_histogram(aes(rbegin1-re.closest1),binwidth=1)+
    scale_x_continuous(breaks=-window:window)
  pright=ggplot(data[abs(rbegin2-re.closest2)<=window])+geom_histogram(aes(rbegin2-re.closest2),binwidth=1)+
    scale_x_continuous(breaks=-window:window)
  pdiag=ggplot(data[abs(rbegin2-rbegin1)<maxlen&strand1==1&strand2==0])+geom_histogram(aes(rbegin2-rbegin1),binwidth=10)
  return(list(pleft=pleft,pright=pright,pdiag=pdiag))
}


#' Wrapper function for the whole preprocessing
#' 
#' See \code{\link{read_tsv}}, \code{\link{categorize_by_new_type}} and
#' \code{\link{prepare_for_sparse_cs_norm}} for further detail and arguments. 
#' 
#' @param infile character. Path to TadBit .tsv file
#' @param outprefix character. Prefix to output intermediate files.
#' @param condition character. WT, KO etc.
#' @param replicate character. Replicate number
#' @param name character. A name for the experiment. By default, it is condition, enzyme and replicate
#' @param enzyme character. HindIII, DpnII etc.
#' @param experiment character. Hi-C Capture-C etc. For now it is Hi-C and cannot be changed
#' @inheritParams read_tsv
#' @inheritParams categorize_by_new_type
#' @inheritParams prepare_for_sparse_cs_norm
#' @param save.data boolean. Whether to save a CSdata object that contains the raw data
#'
#' @return A CSdata object
#' @export
#' 
#' @examples
read_and_prepare = function(infile, outprefix, condition, replicate, enzyme = "HindIII", experiment = "Hi-C",
                            name = paste(condition, enzyme, replicate), skip = 0L, nrows = -1L,
                            circularize = -1, dangling.L = c(0, 4), dangling.R = c(3, -1), maxlen = 600,
                            save.data=T) {
  match.arg(experiment)
  message("*** READ")
  data=read_tsv(infile, skip=skip, nrows=nrows)
  message("*** CATEGORIZE")
  data = categorize_by_new_type(data, dangling.L = dangling.L, dangling.R = dangling.R, maxlen = maxlen)
  message("*** BIASES AND COUNTS")
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=circularize)
  dset_statistics(cs_data$biases,cs_data$counts)
  message("*** WRITE")
  csd = new("CSdata", info=list(name=name, condition=condition, replicate=replicate,
                               enzyme=enzyme, experiment=experiment,
                               dangling.L=deparse(dangling.L), dangling.R=deparse(dangling.R), maxlen=maxlen,
                               filename=infile),
                     settings=list(circularize=circularize),
                     data=data, biases=cs_data$biases, counts=cs_data$counts)
  if (save.data==T) save(csd, file=paste0(outprefix,"_csdata_with_data.RData"))
  csd@data=data.table()
  save(csd, file=paste0(outprefix,"_csdata.RData"))
  return(csd)
}

#' Merge one or more CSdata objects into a CSnorm object
#'
#' @param datasets list of CSdata objects
#'
#' @return CSnorm object
#' @export
#'
#' @examples
merge_cs_norm_datasets = function(datasets) {
  #compile table of experiments, sorted by id
  experiments = rbindlist(lapply(datasets, function(x) x@info))
  setkey(experiments, name)
  stopifnot(experiments[,.N]==experiments[,uniqueN(name)]) #id must be unique
  #all experiments must have the same settings
  stopifnot(uniqueN(rbindlist(lapply(datasets, function(x) x@settings)))==1)
  #merge biases and counts into data tables, make IDs unique
  biases = lapply(datasets, function(x) x@biases)
  counts = lapply(datasets, function(x) x@counts)
  names(biases) <- sapply(datasets, function(x) x@info$name)
  names(counts) <- sapply(datasets, function(x) x@info$name)
  stopifnot(all(sapply(datasets, function(x) x@biases[,all(id==1:.N)])))
  sizes=sapply(biases, function(x) x[,.N])
  sizes=cumsum(sizes)
  if (length(datasets)>1) {
    for (i in 2:length(datasets)) {
      biases[[i]][,id:=id+sizes[i-1]]
      counts[[i]][,id1:=id1+sizes[i-1]]
      counts[[i]][,id2:=id2+sizes[i-1]]
    }
  }
  biases=rbindlist(biases, use.names=T, idcol="name")
  counts=rbindlist(counts, use.names=T, idcol="name")
  stopifnot(biases[,all(id==1:.N)])
  #return CSnorm object
  new("CSnorm", experiments=experiments,
                design=data.table(),
                settings=datasets[[1]]@settings,
                biases=biases, counts=counts)
}
