#' @include csnorm.R
NULL

#' Read TADBit tsv file and return (paired-end) reads as data.table
#'
#' @param fname The filename
#' @param nrows,skip see \code{\link[data.table]{fread}}
#' @param locus if not NULL, a vector in the format c("chr3",123,456) to keep only
#' paired-end reads where both ends map in region chr3:123-456. Otherwise, read everything.
#'
#' @return a data.table
#' @export
#'
#' @examples
read_tsv = function(fname, nrows=-1L, skip=0L, locus=NULL) {
  data = fread(fname, col.names=c("id", "chr1","begin1","strand1", "length1",
                                  "re.up1","re.dn1",
                                  "chr2","begin2","strand2", "length2",
                                  "re.up2","re.dn2"), nrows=nrows, skip=skip)
  if (!is.null(locus)) data=data[chr1==locus[1]&chr2==locus[1]
                                 &begin1>=as.integer(locus[2])&begin2>=as.integer(locus[2])
                                 &begin1<=as.integer(locus[3])&begin2<=as.integer(locus[3])]
  #only one chromosome for now
  stopifnot(data[,nlevels(factor(chr1))]==1)
  stopifnot(data[,nlevels(factor(chr2))]==1)
  stopifnot(data[,chr1[1]==chr2[1]])
  data[,chr1:=NULL]
  data[,chr2:=NULL]
  #
  cat("put read2 always downstream of read1\n")
  bdata=data[begin1>begin2]
  setnames(bdata, c("begin1","strand1", "length1", "re.up1","re.dn1",
                    "begin2","strand2", "length2", "re.up2","re.dn2"),
           c("begin2","strand2", "length2", "re.up2","re.dn2",
             "begin1","strand1", "length1", "re.up1","re.dn1"))
  data = rbind(data[begin1<=begin2],bdata)
  #
  cat("set rbegin and rend\n")
  data[strand1==0,c("rbegin1","rend1"):=list(begin1,begin1-length1+1)]
  data[strand1==1,c("rbegin1","rend1"):=list(begin1,begin1+length1-1)]
  data[strand2==0,c("rbegin2","rend2"):=list(begin2,begin2-length2+1)]
  data[strand2==1,c("rbegin2","rend2"):=list(begin2,begin2+length2-1)]
  #
  cat("find each read's closest restriction site\n")
  data[,re.closest1:=ifelse(rbegin1-re.up1 < re.dn1-rbegin1,re.up1,re.dn1)]
  data[,re.closest2:=ifelse(rbegin2-re.up2 < re.dn2-rbegin2,re.up2,re.dn2)]
  #
  cat("number it\n")
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
  cat("sort by begins\n")
  setkey(data,rbegin1,rbegin2,strand1,strand2)
  #
  cat("remove duplicates\n")
  data=unique(data)
  return(data)
}

#' Get subset of raw reads
#'
#' @param data reads data.table containing the raw reads
#' @param b1,e1 begin / end of read range you want
#' @param b2,e2 (optional) if provided, give extradiagonal portion 
#'
#' @return reads data.table subset
#' @export
#'
#' @examples
get_raw_reads = function(data, b1, e1, b2=NULL, e2=NULL) {
  if (is.null(b2)) b2=b1
  if (is.null(e2)) e2=e1
  stopifnot(e1>b1)
  stopifnot(e2>b2)
  stopifnot(b2>=b1)
  return(data[begin1>=b1 & begin1+length1<=e1 & begin2>=b2 & begin2+length2<=e2])
}

#' Bin a reads data.table into a matrix of a given resolution
#'
#' @return a data.table representing the binned data
#' @keywords internal
#'
#' @examples

bin_data = function(data, resolution, b1=NULL, b2=NULL) {
  if (is.null(b1)) b1=data[,min(begin1)]-1
  if (is.null(b2)) b2=data[,min(begin2)]-1
  bins1=seq(b1,data[,max(begin1)]+resolution,resolution)
  bins2=seq(b2,data[,max(begin2)]+resolution,resolution)
  #
  sub = data[,.(begin1,begin2,bin1=cut2(begin1, bins1, oneval=F, onlycuts=T, digits=10),
                bin2=cut2(begin2, bins2, oneval=F, onlycuts=T, digits=10))
             ][,.N,by=c("bin1","bin2")]
  #
  sub[,begin1:=do.call(as.integer, tstrsplit(as.character(bin1), "[[,]")[2])]
  sub[,end1:=do.call(as.integer, tstrsplit(as.character(bin1), "[],)]")[2])]
  sub[,begin2:=do.call(as.integer, tstrsplit(as.character(bin2), "[[,]")[2])]
  sub[,end2:=do.call(as.integer, tstrsplit(as.character(bin2), "[],)]")[2])]
  return(sub)
}

#' Plot reads in binned form at a given resolution
#'
#' @param dt a reads data.table containing the parsed tsv reads
#' @param resolution the requested resolution, in base pairs
#' @param b1,e1 the begin and end of the range requested 
#' @param b2,e2 optional, if you want an outer-diagonal section 
#' @param diagonal logical, default FALSE. Should the diagonal be indicated? 
#' @param rsites logical, default FALSE. Should restriction sites be indicated?
#'
#' @return
#' @export
#'
#' @examples
plot_binned = function(dt, resolution, b1, e1, b2=NULL, e2=NULL, diagonal=F, rsites=F) {
  if (is.null(b2)) b2=b1
  if (is.null(e2)) e2=e1
  binned = bin_data(dt, resolution, b1=b1, b2=b2)
  #
  p=ggplot(binned, aes(begin1,begin2, fill=log(N)))+geom_raster()+
    coord_cartesian(xlim=c(b1,e1),ylim=c(b2,e2))+
    scale_fill_gradient(low="white", high="black")
  if (diagonal) p = p+stat_function(fun=function(x){x})
  if (rsites) {
    rsites1=dt[,unique(c(re.up1,re.dn2))]
    rsites1=rsites1[rsites1>=b1&rsites1<=e1]
    rsites2=dt[,unique(c(re.up2,re.dn2))]
    rsites2=rsites2[rsites2>=b1&rsites2<=e1]
    p = p + geom_vline(xintercept=rsites1) + geom_hline(yintercept=rsites2)
  }
  print(p)
}

#' Plot raw reads
#'
#' @inheritParams plot_binned
#'
#' @return
#' @export
#'
#' @examples
plot_raw = function(dt, b1=NULL, e1=NULL, b2=NULL, e2=NULL, diagonal=T, rsites=T) {
  if (is.null(b1)) b1 = dt[,min(rbegin1)]
  if (is.null(e1)) e1 = dt[,max(rend2)]
  faceted = (!is.null(b2)) & (!is.null(e2)) & b2>=e1 & diagonal==T #show 3 viewpoints simultaneously
  if (length(faceted)!=1) faceted=F #happens if b2 is NULL
  #get subset of data, add facet info
  if (faceted==F) {
    sub = get_raw_reads(dt, b1, e1, b2, e2)
    if (is.null(b2)) b2=b1
    if (is.null(e2)) e2=e1
    p=ggplot(sub[,.SD[sample(.N,min(.N,100000))],,by=c("strand1","strand2")])
    if (diagonal==T) p = p+stat_function(fun=function(x){x})
    if (rsites==T) {
      #rsite positions
      rsites.all=sub[,unique(c(re.up1,re.up2,re.dn1,re.dn2))]
      rsites1=rsites.all[rsites.all>=b1&rsites.all<=e1]
      p = p + geom_vline(xintercept=rsites1, colour="grey")
      rsites2=rsites.all[rsites.all>=b2&rsites.all<=e2]
      p = p + geom_hline(yintercept=rsites2, colour="grey")
    }
    #plot
    p=p+geom_segment(aes(x=rbegin1, y=rbegin2, xend=rend1, yend=rend2, colour=category),
                     arrow=arrow(type="closed", length=unit(0.005,"npc")))
    p=p+xlim(b1,e1)+ylim(b2,e2)
  } else {
    sub1 = get_raw_reads(dt, b1, e1)
    sub1[,c("xfac","yfac"):=list(1,0)]
    sub2 = get_raw_reads(dt, b2, e2)
    sub2[,c("xfac","yfac"):=list(0,1)]
    sub3 = get_raw_reads(dt, b1, e1, b2, e2)
    sub3[,c("xfac","yfac"):=list(0,0)]
    sub=rbind(sub1,sub2,sub3)
    #diagonal
    diag.data = data.table(x=c(b1,b2), y=c(b1,b2), xend=c(e1,e2), yend=c(e1,e2), xfac=c(1,0), yfac=c(0,1))
    p=ggplot(sub[,.SD[sample(.N,min(.N,100000))],,by=c("strand1","strand2")])
    p = p+geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=diag.data)
    #rsites
    if (rsites==T) {
      rsites.all=sub[,unique(c(re.up1,re.up2,re.dn1,re.dn2))]
      rsites1=rsites.all[rsites.all>=b1&rsites.all<=e1]
      rsites2=rsites.all[rsites.all>=b2&rsites.all<=e2]
      rsite.data.x = data.table(z=c(rsites1,rsites1,rsites2),
                                xfac=c(rep(c(1,0),each=length(rsites1)), rep(0,length(rsites2))),
                                yfac=c(rep(0,length(rsites1)*2), rep(1,length(rsites2))))
      rsite.data.y = data.table(z=c(rsites2,rsites2,rsites1),
                                yfac=c(rep(c(1,0),each=length(rsites2)), rep(0,length(rsites1))),
                                xfac=c(rep(0,length(rsites2)*2), rep(1,length(rsites1))))
      p = p + geom_vline(aes(xintercept=z), rsite.data.x, colour="grey") + geom_hline(aes(yintercept=z), rsite.data.y, colour="grey")
    }
    #plot
    p=p+geom_segment(aes(x=rbegin1, y=rbegin2, xend=rend1, yend=rend2, colour=category),
                     arrow=arrow(type="closed", length=unit(0.005,"npc")))
    p=p+facet_grid(xfac~yfac, scales = "free", space="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  print(p)
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
#' @param read.len integer vector. The length of a full read, used for dangling ends (temporary)
#'   
#' @return The same data.table, with an additional "category" column
#' @keywords internal
#' 
#' @examples
categorize_by_new_type = function(sub, maxlen=600, dangling.L = c(0), dangling.R = c(3), read.len=40) {
  sub[,category:="other"]
  #careful: order is important because attribution criteria overlap
  cat("Contacts\n")
  sub[re.closest1.idx != re.closest2.idx & re.dn1 != re.dn2,
      category:=ifelse(strand1==1 & re.dn1 - rbegin1 < maxlen & strand2==1 & re.dn2 - rbegin2 < maxlen, "contact up",
                       ifelse(strand1==1 & re.dn1 - rbegin1 < maxlen & strand2==0 & rbegin2 - re.up2 < maxlen & rbegin2-rbegin1 >= maxlen, "contact far", #necessary to not overwrite DE
                              ifelse(strand1==0 & rbegin1 - re.up1 < maxlen & strand2==1 & re.dn2 - rbegin2 < maxlen, "contact close",
                                     ifelse(strand1==0 & rbegin1 - re.up1 < maxlen & strand2==0 & rbegin2 - re.up2 < maxlen, "contact down", "other"))))]
  cat("Random\n")
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & re.dn1 == re.dn2, category:="random"]
  cat("Rejoined\n")
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & re.dn1 <= re.up2, category:="rejoined"]
  sub[category=="rejoined" & ( !(length1 %in% read.len) | !(length2 %in% read.len) ), category:="random"]
  cat("Self Circle\n")
  sub[re.dn1 == re.dn2 & strand1==0 & strand2==1 & rbegin1 - re.up1 < maxlen & re.dn2 - rbegin2 < maxlen, category:="self circle"]
  cat("Dangling L/R\n")
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & ((rbegin1 - re.closest1) %in% dangling.L) & length1 %in% read.len, category:="dangling L"]
  sub[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & ((rbegin2 - re.closest2) %in% dangling.R) & length2 %in% read.len, category:="dangling R"]
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
#' @keywords internal
#' 
#' @examples
prepare_for_sparse_cs_norm = function(data, both=F, circularize=-1) {
  #get counts that are not other, random or self circle
  cat("Sum up counts and biases\n")
  enrich = data[!(category %in% c("other","random", "self circle")), .N,
                keyby=c("re.closest1","re.closest2", "category")]
  enrich[,category:=gsub(" ", ".", category)]
  #
  #dangling ends: attribute to cut site
  cat("Reattribute biases to proper rsites\n")
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
  #cat("RSites list\n")
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
  cat("Add Rsite info\n")
  enrich = merge(enrich, rsites, by.x="re.closest1", by.y="re.pos")
  enrich = merge(enrich, rsites, by.x="re.closest2", by.y="re.pos", suffixes=c("1", "2"))
  #
  cat("Produce bias list\n")
  X=dcast.data.table(enrich[id1 == id2 & category %in% c("dangling.L", "dangling.R", "rejoined"),
                            .(id=id1, pos=re.closest1, category, N)],
                     ...~category, value.var="N", fun.aggregate = sum)
  X=rbind(X,rsites[!(id%in%X[,id]),.(id,pos=re.pos,dangling.L=0,dangling.R=0,rejoined=0)])
  setkey(X,id)
  stopifnot(all(X[,id==.I])) #IDs must start at one and have no gaps, for current stan implementation
  #
  cat("Produce count list\n")
  Y=dcast.data.table(enrich[id1 != id2 
                            & category %in% c("contact.close", "contact.down", "contact.far", "contact.up"),
                            .(id1,id2,pos1=re.closest1,pos2=re.closest2,category,N)],
                     ...~category, value.var="N", fill=0)
  setkey(Y,id1,id2)
  Y[,distance:=pos2-pos1]
  if (circularize>0) Y[,distance:=pmin(distance,circularize-distance+1)]
  #
  if (both == T) {
    cat("Merge them\n")
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
#' @keywords internal
#'
dset_statistics = function(biases,counts){
  cat("Mean counts\n")
  cat("   Rejoined  : ", biases[,mean(rejoined)],"\n")
  cat("   Dangling L: ", biases[,mean(dangling.L)],"\n")
  cat("   Dangling R: ", biases[,mean(dangling.R)],"\n")
  cat("   C. close  : ", counts[,mean(contact.close)],"\n")
  cat("   C. far    : ", counts[,mean(contact.far)],"\n")
  cat("   C. up     : ", counts[,mean(contact.up)],"\n")
  cat("   C. down   : ", counts[,mean(contact.down)],"\n")
  cat("Median counts\n")
  cat("   Rejoined  : ", biases[,median(rejoined)],"\n")
  cat("   Dangling L: ", biases[,median(dangling.L)],"\n")
  cat("   Dangling R: ", biases[,median(dangling.R)],"\n")
  cat("   C. close  : ", counts[,median(contact.close)],"\n")
  cat("   C. far    : ", counts[,median(contact.far)],"\n")
  cat("   C. up     : ", counts[,median(contact.up)],"\n")
  cat("   C. down   : ", counts[,median(contact.down)],"\n")
  cat("Percent of zero counts\n")
  cat("   Rejoined  : ", biases[rejoined==0,.N]/biases[,.N]*100,"\n")
  cat("   Dangling L: ", biases[dangling.L==0,.N]/biases[,.N]*100,"\n")
  cat("   Dangling R: ", biases[dangling.R==0,.N]/biases[,.N]*100,"\n")
  cat("   C. close  : ", counts[contact.close==0,.N]/counts[,.N]*100,"\n")
  cat("   C. far    : ", counts[contact.far==0,.N]/counts[,.N]*100,"\n")
  cat("   C. up     : ", counts[contact.up==0,.N]/counts[,.N]*100,"\n")
  cat("   C. down   : ", counts[contact.down==0,.N]/counts[,.N]*100,"\n")
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
generate_fake_dataset = function(biases.ref=NULL, num_rsites=3000, genome_size=1000000, eC=-4.4, eRJ=3, eDE=1, alpha=2,
                                 condition="WT", replicate="1", enzyme="NA",
                                 name = paste("Fake", condition, enzyme, replicate), dmin=1000, signal=F) {
  if (is.null(biases.ref)) {
    #place rsites
    biases=data.table(id=seq(num_rsites),
                      pos=cumsum(1+rmultinom(n=1, size=genome_size-num_rsites,
                                                  prob=rep(1,num_rsites+1)))[1:num_rsites])
    setkey(biases,id)
  } else {
    biases=biases.ref[,.(id,pos=(pos-min(pos)+1))]
    num_rsites=biases[,.N]
    genome_size=biases[,max(pos)+1]
  }
  #build biases
  biases[,true_log_iota:=(pos-genome_size/2)/genome_size]
  biases[,true_log_iota:=true_log_iota-mean(true_log_iota)]
  #ggplot(biases,aes(pos,true_log_iota))+geom_point()+geom_line()
  biases[,true_log_rho:=sin(pos/1300)]
  biases[,true_log_rho:=true_log_rho-mean(true_log_rho)]
  #ggplot(biases,aes(pos,true_log_rho))+geom_point()+geom_line()
  biases[,true_log_mean_RJ:=eRJ+(true_log_iota+true_log_rho)/2]
  biases[,true_log_mean_DL:=eDE+true_log_iota]
  biases[,true_log_mean_DR:=eDE+true_log_rho]
  #draw dangling/rejoined
  biases[,dangling.L:=rnbinom(.N, mu=exp(true_log_mean_DL), size=alpha)]
  biases[,dangling.R:=rnbinom(.N, mu=exp(true_log_mean_DR), size=alpha)]
  biases[,rejoined:=rnbinom(.N, mu=exp(true_log_mean_RJ), size=alpha)]
  #report rsites in counts
  counts=CJ(id1=biases[,id],id2=biases[,id])[id2>id1]
  counts=merge(counts, biases[,.(id1=id,pos,true_log_iota,true_log_rho)], by="id1")
  counts=merge(counts, biases[,.(id2=id,pos,true_log_iota,true_log_rho)], by="id2", suffixes=c("1","2"))
  setkey(counts, id1, id2)
  #build decay
  counts[,distance:=pos2-pos1]
  counts[,true_log_decay:=100*dnorm(log(distance), mean=0, sd=log(max(distance))/3)]
  counts[,true_log_decay:=true_log_decay-mean(true_log_decay)]
  #ggplot(counts[sample(.N,min(100000,.N))])+geom_line(aes(distance,true_log_decay))+scale_x_log10()
  counts[,base_count:=eC+true_log_decay]
  #add signal
  if (signal==T) {
    #loop
    loop_center=c(genome_size*0.2, genome_size*0.8)
    loop_sd=c(genome_size/50,genome_size/20)
    counts[,true_phi:=2*exp(-(pos1-loop_center[1])^2/loop_sd[1]^2)*exp(-(pos2-loop_center[2])^2/loop_sd[2]^2)]
    #tad
    tad_begin=genome_size*0.6
    tad_end=genome_size*0.8
    counts[,true_phi:=true_phi+ifelse(pos1>=tad_begin&pos1<=tad_end&pos2>=tad_begin&pos2<=tad_end,1,0)]
  } else {
    counts[,true_phi:=0]
  }
  counts[,true_log_mean_cclose:=base_count+true_phi+true_log_rho1+true_log_iota2]
  counts[,true_log_mean_cfar:=base_count+true_phi+true_log_iota1+true_log_rho2]
  counts[,true_log_mean_cup:=base_count+true_phi+true_log_iota1+true_log_iota2]
  counts[,true_log_mean_cdown:=base_count+true_phi+true_log_rho1+true_log_rho2]
  #draw counts
  counts[,contact.close:=rnbinom(.N, mu=exp(true_log_mean_cclose), size=alpha)]
  counts[,contact.far:=rnbinom(.N, mu=exp(true_log_mean_cfar), size=alpha)]
  counts[,contact.up:=rnbinom(.N, mu=exp(true_log_mean_cup), size=alpha)]
  counts[,contact.down:=rnbinom(.N, mu=exp(true_log_mean_cdown), size=alpha)]
  #
  #counts[,bin1:=round(pos1/10000)]
  #counts[,bin2:=round(pos2/10000)]
  #binned=counts[,.(count=sum(contact.far+contact.close+contact.up+contact.down)),by=c("bin1","bin2")]
  #ggplot(binned)+geom_raster(aes(bin1,bin2,fill=log(count)))
  #statistics
  #csnorm:::dset_statistics(biases,counts)
  counts=counts[contact.close+contact.far+contact.up+contact.down>0]
  #csd object
  dmax=biases[,max(pos)-min(pos)]+0.01
  csd = new("CSdata", info=list(name=name, condition=condition, replicate=replicate,
                                enzyme=enzyme, experiment="Hi-C",
                                dangling.L="NA", dangling.R="NA", maxlen="NA",
                                filename="NA"),
            settings=list(circularize=-1, dmin=dmin, dmax=dmax),
            data=data.table(), biases=biases, counts=counts)
  return(csd)
}

#' Diagnostic plots to determine dangling.L, dangling.R and maxlen arguments
#' 
#' @param infile A tadbit .tsv file as character string, or a data.table returned by read.tsv
#' @param window how many bases to plot around cut site
#' @param maxlen how many bases to plot off-diagonal
#' @param skip.fbm boolean. If TRUE (default), skip fragment-based IDs, 
#'   containing # and ~. Otherwise keep them.
#' @param read.len length of reads admissible for dangling end candidates
#' @inheritParams read_tsv
#'   
#' @return Two plots, pdangling to determine dangling.L and dangling.R, and
#'   pdiag for maxlen. Also returns data, the reads data table used for the
#'   plots.
#' @export
#' 
#' @examples
examine_dataset = function(infile, window=15, maxlen=1000, skip.fbm=T, read.len=40, skip=0L, nrows=-1L, locus=NULL) {
  if (is.character(infile)) {
    data=read_tsv(infile, skip=skip, nrows=nrows, locus=locus)
  } else {
    data=infile
  }
  if (skip.fbm == T) data = data[!grepl("[#~]",id)]
  #dangling
  dleft=data[abs(rbegin1-re.closest1)<=window&abs(rbegin2-rbegin1)<maxlen&strand1==1&strand2==0&length1%in%read.len,
             .(dist=rbegin1-re.closest1,cat="left")]
  dright=data[abs(rbegin2-re.closest2)<=window&abs(rbegin2-rbegin1)<maxlen&strand1==1&strand2==0&length2%in%read.len,
              .(dist=rbegin2-re.closest2,cat="right")]
  pdangling=ggplot(rbind(dleft,dright))+
    geom_histogram(aes(dist),binwidth=1)+scale_x_continuous(breaks=-window:window)+facet_grid(~cat)+
    xlab("distance from cut site")+ylab("number of reads")
  #read size distribution
  pdiag=ggplot(data[abs(rbegin2-rbegin1)<maxlen&strand1==1&strand2==0])+
    geom_histogram(aes(rbegin2-rbegin1),binwidth=10)+
    xlab("sonication fragment size")+ylab("number of reads")
  #close decay
  dclose=data[abs(re.closest1-re.closest2)<3000 & re.closest1!=re.closest2 & re.dn1 != re.dn2]
  dclose=dclose[,.(count=.N,dist=abs(re.closest2-re.closest1)[1]),by=c("re.closest1.idx","re.closest2.idx")]
  dclose[,dbin:=cut(dist,1000,ordered_result=T,right=F,include.lowest=T,dig.lab=12)]
  pclose=ggplot(dclose[,.(dist=mean(dist),N=mean(count)),by=dbin])+geom_point(aes(dist,N))+
    xlab("mean distance between cut sites")+ylab("mean count between cut sites")
  return(list(pdangling=pdangling,pdiag=pdiag,pclose=pclose,data=data))
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
#' @param name character. A name for the experiment. By default, it is
#'   condition, enzyme and replicate
#' @param enzyme character. HindIII, DpnII etc.
#' @param experiment character. Hi-C Capture-C etc. For now it is Hi-C and
#'   cannot be changed
#' @inheritParams read_tsv
#' @inheritParams categorize_by_new_type
#' @inheritParams prepare_for_sparse_cs_norm
#' @param dmin numeric. Minimum disrance between cut sites to consider this as a
#'   contact. Note that this setting will only be enforced upon merging, and the
#'   returned object might contain contacts at a smaller distance.
#' @param save.data boolean. Whether to save a CSdata object that contains the
#'   raw data
#'   
#' @return A CSdata object
#' @export
#' 
#' @examples
read_and_prepare = function(infile, outprefix, condition, replicate, enzyme = "HindIII", experiment = "Hi-C",
                            name = paste(condition, enzyme, replicate), skip = 0L, nrows = -1L, locus=NULL,
                            circularize = -1, dangling.L = c(0, 4), dangling.R = c(3, -1), maxlen = 600,
                            read.len=40, dmin=2000, save.data=T) {
  match.arg(experiment)
  cat("*** READ\n")
  data=read_tsv(infile, skip=skip, nrows=nrows, locus=locus)
  if (circularize>0) data=data[!((re.closest1 %in% c(1,circularize)) | (re.closest2 %in% c(1,circularize)))]
  cat("*** CATEGORIZE\n")
  data = categorize_by_new_type(data, dangling.L = dangling.L, dangling.R = dangling.R, maxlen = maxlen, read.len=read.len)
  cat("*** BIASES AND COUNTS\n")
  cs_data = prepare_for_sparse_cs_norm(data, both=F, circularize=circularize)
  dset_statistics(cs_data$biases,cs_data$counts)
  cat("*** WRITE\n")
  #determine dmax
  if (circularize>0) {
    dmax=circularize/2+0.01
  } else {
    dmax=cs_data$biases[,max(pos)-min(pos)]+0.01
  }
  csd = new("CSdata", info=list(name=name, condition=condition, replicate=replicate,
                                enzyme=enzyme, experiment=experiment,
                                dangling.L=deparse(dangling.L), dangling.R=deparse(dangling.R), maxlen=maxlen,
                                filename=infile),
            settings=list(circularize=circularize, dmin=dmin, dmax=dmax),
            data=data, biases=cs_data$biases, counts=cs_data$counts)
  if (save.data==T) save(csd, file=paste0(outprefix,"_csdata_with_data.RData"))
  csd@data=data.table()
  save(csd, file=paste0(outprefix,"_csdata.RData"))
  return(csd)
}

#' Fuse cut sites that are very close to each other
#'
#' @param csd CSData object
#' @keywords internal
#' @return a list with filtered biases and counts
#'
fuse_close_cut_sites = function(biases,counts,dfuse,name,circularize) {
  #group cut sites
  fused=biases[,.(id,pos,prev.pos=pos-shift(pos,type = "lag"),incr.id=ifelse(pos-shift(pos,type = "lag")>dfuse,1,0))]
  fused[1,incr.id:=1]
  fused[,newid:=cumsum(incr.id)]
  setkey(fused,id)
  fused[,newpos:=as.integer(mean(pos)),by=newid]
  fused[,c("prev.pos","incr.id"):=NULL]
  cat(name, ": fused ", fused[,.N]-fused[,uniqueN(newid)], " cut sites (",
      (1-fused[,uniqueN(newid)]/fused[,.N])*100, " %)\n")
  #transfer new IDs
  biases=fused[biases,.(id=newid,pos=newpos,dangling.L,dangling.R,rejoined)]
  counts=merge(counts[,.(id1,id2,contact.close,contact.down,contact.far,contact.up)],fused,by.x="id1",by.y="id")
  counts=merge(counts,fused,by.x="id2",by.y="id",suffixes=c("1","2"))
  counts[,c("id1","id2","pos1","pos2"):=NULL]
  setnames(counts,c("newid1","newid2","newpos1","newpos2"),c("id1","id2","pos1","pos2"))
  #sum counts and average positions
  biases=biases[,.(pos=pos[1],dangling.L=sum(dangling.L),
                   dangling.R=sum(dangling.R),rejoined=sum(rejoined)),keyby=id]
  counts=counts[id1!=id2,.(pos1=pos1[1],pos2=pos2[1],
                   contact.close=sum(contact.close),contact.down=sum(contact.down),
                   contact.far=sum(contact.far),contact.up=sum(contact.down)),
                keyby=c("id1","id2")]
  counts[,distance:=pos2-pos1]
  if (circularize>0) counts[,distance:=pmin(distance,circularize-distance+1)]
  return(list(biases=biases,counts=counts))
}

#' Merge one or more CSdata objects into a CSnorm object and set experimental design
#'
#' This will do some cleanup on the data, namely:
#' - remove all contacts lower than their respective dmin
#' - remove the counts at cut sites with too low or too high coverage (default 1st and 99th percentile)
#' Note that you can customize the design, as long as you keep the numbering sequential
#' and starting at 1.
#' 
#' @param datasets list of CSdata objects
#' @param different.decays character. If "none" (default), one decay is modelled
#'   for all experiments. If "all", each experiment has its own decay. If
#'   "enzyme", experiments with a different enzyme get their own decay. If
#'   "condition", experiments with a different condition get their own decay.
#'   Arguments can be abbreviated and combined
#' @param dfuse Fuse cut sites that are closer than dfuse to each other
#'
#' @return CSnorm object
#' @export
#'
#' @examples
merge_cs_norm_datasets = function(datasets, different.decays=c("none","all","enzyme","condition"), dfuse=5, qmin=0.01) {
  cat("compile table of experiments, sorted by id\n")
  experiments = rbindlist(lapply(datasets, function(x) x@info))
  experiments[,name:=ordered(name, levels=name)]
  setkey(experiments, name)
  if (!(experiments[,.N]==experiments[,uniqueN(name)])) stop("experiment names must be unique")
  if (uniqueN(rbindlist(lapply(datasets, function(x) x@settings[c("dmin","circularize")])))!=1)
    stop("all experiments must have the same dmin and circularize settings!")
  settings=datasets[[1]]@settings[c("dmin","circularize")]
  settings$dmax=max(sapply(datasets, function(x) x@settings$dmax))
  settings$dfuse=dfuse
  settings$qmin=qmin
  #
  filtered = lapply(datasets, function(csd){fuse_close_cut_sites(csd@biases, csd@counts, dfuse,
                                                                 csd@info$name, csd@settings$circularize)})
  biases = lapply(filtered, function(x) x$biases)
  counts = lapply(filtered, function(x) x$counts)
  names(biases) <- sapply(datasets, function(x) x@info$name)
  names(counts) <- sapply(datasets, function(x) x@info$name)
  stopifnot(all(sapply(datasets, function(x) x@biases[,all(id==1:.N)])))
  #
  cat("merge biases and counts into data tables, make IDs unique\n")
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
  biases[,name:=ordered(name,levels=experiments[,levels(name)])]
  setkey(biases,name,id,pos)
  counts[,name:=ordered(name,levels=experiments[,levels(name)])]
  setkey(counts,name,id1,pos1,id2,pos2)
  #
  cat("design matrix\n")
  different.decays=match.arg(different.decays, several.ok = T)
  design=experiments[,.(name,enzyme,condition)]
  design[,genomic:=as.integer(factor(enzyme))]
  if (setequal(different.decays,"none")) {
    design[,decay:=1]
  } else if (setequal(different.decays,"all")) {
    design[,decay:=.I]
  } else if (setequal(different.decays,"enzyme")) {
    design[,decay:=as.integer(factor(enzyme))]
  } else if (setequal(different.decays,"condition")) {
    design[,decay:=as.integer(factor(condition))]
  } else if (setequal(different.decays, c("condition","enzyme"))) {
    design[,decay:=as.integer(factor(paste(condition,enzyme)))]
  } else {
    stop("I do not know what to do with different.decays = ", different.decays)
  }
  design[,c("enzyme","condition"):=list(NULL,NULL)]
  cat("enforce minimum distance\n")
  counts=counts[distance>=settings$dmin]
  cat("return CSnorm object\n")
  new("CSnorm", experiments=experiments,
      design=design, zeros=data.table(),
      settings=settings,
      biases=biases, counts=counts)
}

#' Take a random subset of the reads of a given csnorm object
#' 
#' @param csnorm a csnorm object
#' @param subsampling.pc positive integer
#'   
#' @return The csnorm object with updated biases and counts slots
#' @export
#' 
#' @examples
subsample_csnorm = function(cs, subsampling.pc=100) {
  biases=copy(cs@biases)
  biases[,dangling.L:=rbinom(.N,dangling.L,subsampling.pc/100)]
  biases[,dangling.R:=rbinom(.N,dangling.R,subsampling.pc/100)]
  biases[,rejoined:=rbinom(.N,rejoined,subsampling.pc/100)]
  counts=copy(cs@counts)
  counts[,contact.close:=rbinom(.N,contact.close,subsampling.pc/100)]
  counts[,contact.down:=rbinom(.N,contact.down,subsampling.pc/100)]
  counts[,contact.far:=rbinom(.N,contact.far,subsampling.pc/100)]
  counts[,contact.up:=rbinom(.N,contact.up,subsampling.pc/100)]
  new("CSnorm", experiments=cs@experiments,
      design=cs@design,
      settings=cs@settings,
      biases=biases, counts=counts)
}

