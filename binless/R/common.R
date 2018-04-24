#' @include binless.R
NULL

#' Compute means for positive and zero counts using previous params
#' 
#' @return data.table containing the following columns
#' - name: name of the dataset
#' - id1, pos1, bin1: coordinates of the cut site
#' - bin2, dbin: signal/distance bin in which we looked at the intersections
#' - dir: fwd (rev) for contacts with downstream (resp upstream) cut-sites
#' - cat: whether we consider contacts on the left (contact L) or on the right (contact R) of this cut site
#' - count: the value of the count in this signal/distance/direction/category bin.
#'  We discard anything below cs@settings$dmin. Note that sum(count) is twice the total number of counts (per dataset).
#' - lmu.nosig: the log(mean) of the background model (exposure+bias+decay)
#' - phi: the log(signal)
#' - mu: the mean, including signal
#' - nobs: how many observations in this signal/distance/direction/category bin have this count.
#'  Note that sum(nobs) is twice the total number of observable counts (per dataset).
#' - eC, log_decay, log_bias: exposure, log(decay) and log(bias) of this dataset. The bias is either iota or rho depending on cat.
#' - z: count/mu-1
#' - var: 1/mu+1/dispersion
#'   
#' @keywords internal
#' 
gauss_common_muhat_mean = function(cs, zeros, sbins) {
  init=cs@par
  ### positive counts (twice, in both directions)
  #compute means
  bsub=merge(init$biases[cat == "contact L",.(genomic.grp=group,pos,log_iota=eta)],
             init$biases[cat == "contact R",.(genomic.grp=group,pos,log_rho=eta)], by=c("genomic.grp","pos"))
  cpos=merge(cbind(cs@design[,.(name,decay.grp=decay,genomic.grp=genomic)],eC=init$eC), cs@counts, by="name",all.x=F,all.y=T)
  setnames(bsub,c("genomic.grp","pos1","log_iota1","log_rho1"))
  cpos = merge(cpos, bsub, all.x=T, all.y=F, by=c("genomic.grp","pos1"))
  setnames(bsub,c("genomic.grp","pos2","log_iota2","log_rho2"))
  cpos = merge(cpos, bsub, all.x=T, all.y=F, by=c("genomic.grp","pos2"))
  setnames(bsub,c("genomic.grp","pos","log_iota","log_rho"))
  cpos[,genomic.grp:=NULL]
  cpos[,c("bin1","bin2","dbin"):=
         list(cut(pos1, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(pos2, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
              cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cpos=init$decay[,.(decay.grp=group,dbin,log_decay)][cpos,,on=c("decay.grp","dbin")]
  cpos[,log_mu.base:=eC + log_decay]
  cpos[,c("lmu.far","lmu.down","lmu.close","lmu.up"):=list(log_mu.base+log_iota1+log_rho2,
                                                           log_mu.base+log_rho1 +log_rho2,
                                                           log_mu.base+log_rho1 +log_iota2,
                                                           log_mu.base+log_iota1+log_iota2)]
  cpos[,log_mu.base:=NULL]
  cpos=rbind(cpos[contact.close>0,.(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd",
                                    count=contact.close, lmu.nosig=lmu.close, nobs=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.far>0,  .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd",
                                    count=contact.far,   lmu.nosig=lmu.far,   nobs=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.down>0, .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact R", dir="fwd",
                                    count=contact.down,  lmu.nosig=lmu.down,  nobs=1, eC, log_decay, log_bias=log_rho1)],
             cpos[contact.up>0,   .(name, id1=id1, pos1=pos1, bin1=bin1, bin2=bin2, dbin, cat="contact L", dir="fwd",
                                    count=contact.up,    lmu.nosig=lmu.up,    nobs=1, eC, log_decay, log_bias=log_iota1)],
             cpos[contact.far>0,  .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev",
                                    count=contact.far,   lmu.nosig=lmu.far,   nobs=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.close>0,.(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev",
                                    count=contact.close, lmu.nosig=lmu.close, nobs=1, eC, log_decay, log_bias=log_iota2)],
             cpos[contact.down>0, .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact R", dir="rev",
                                    count=contact.down,  lmu.nosig=lmu.down,  nobs=1, eC, log_decay, log_bias=log_rho2)],
             cpos[contact.up>0,   .(name, id1=id2, pos1=pos2, bin1=bin2, bin2=bin1, dbin, cat="contact L", dir="rev",
                                    count=contact.up,    lmu.nosig=lmu.up,    nobs=1, eC, log_decay, log_bias=log_iota2)])
  ### zero counts (twice, in both directions)
  czero = cs@design[,.(name,decay.grp=decay,genomic.grp=genomic)][zeros]
  czero = init$decay[,.(decay.grp=group,dbin,log_decay)][czero,.(genomic.grp,name,dbin,bin1,bin2,cat,dir,id1,pos1,nobs=nzero,log_decay),on=c("decay.grp","dbin")]
  stopifnot(czero[is.na(log_decay),.N]==0)
  setnames(bsub,"pos","pos1")
  czero = bsub[czero,on=c("genomic.grp","pos1")]
  czero[,log_bias:=ifelse(cat=="contact L", log_iota, log_rho)]
  czero[,c("log_iota","log_rho"):=NULL]
  czero = cbind(cs@design[,.(name)],eC=init$eC)[czero,on="name"]
  czero[,lmu.nosig:=eC + log_decay + log_bias] #we don't average over j
  czero[,count:=0]
  czero[,genomic.grp:=NULL]
  cts=rbind(cpos,czero)
  setkeyv(cts,key(zeros))
  ### add signal
  signal = binless:::get_signal_matrix(cs, resolution = sbins[2]-sbins[1], groups=cs@experiments[,.(name,groupname=name)])
  signal=rbind(signal[,.(name,bin1,bin2,phi)],signal[bin1!=bin2,.(name,bin1=bin2,bin2=bin1,phi)])
  cts=signal[cts,,on=c("name","bin1","bin2")]
  cts[,mu:=exp(lmu.nosig+phi)]
  ### finalize
  cts[,c("z","var"):=list(count/mu-1,(1/mu+1/init$alpha))]
  setkey(cts,name,id1,pos1,bin1,bin2,dbin,dir,cat)
  #ggplot(cts[name=="T47D es 60 MboI 1"&cat=="contact R"])+geom_line(aes(dbin,log_decay,colour=count>0,group=count>0))
  return(cts)
}

#' Compute predicted means for dangling and rejoined
#' @keywords internal
#' 
gauss_common_muhat_mean_biases = function(cs) {
    #biases
    init=cs@par
    bsub=copy(cs@biases)
    bsub=merge(cbind(cs@design[,.(name,group=genomic)],eRJ=init$eRJ,eDE=init$eDE), bsub, by="name",all.x=F,all.y=T)
    bts=rbind(bsub[,.(group,cat="rejoined",pos, count=rejoined,expo=eRJ,nobs=1)],
              bsub[,.(group,cat="dangling L",pos, count=dangling.L,expo=eDE,nobs=1)],
              bsub[,.(group,cat="dangling R",pos, count=dangling.R,expo=eDE,nobs=1)])
    setkey(bts,group,cat,pos)
    bts = bts[init$biases[cat%in%c("rejoined","dangling L","dangling R"),.(group,cat,pos,eta)]]
    bts[,mu:=exp(expo+eta)]
    return(bts[,.(group,cat,pos,count,mu,nobs,eta)])
}

#' Add begin and end for a binned matrix
#'
#' @keywords internal
#' @export
#'
#' @examples
add_bin_begin_and_end = function(mat) {
  bin1.begin=mat[,bin1]
  bin1.end=mat[,bin1]
  bin2.begin=mat[,bin2]
  bin2.end=mat[,bin2]
  levels(bin1.begin) <- tstrsplit(as.character(levels(bin1.begin)), "[][,)]")[[2]]
  levels(bin1.end) <- tstrsplit(as.character(levels(bin1.end)), "[][,)]")[[3]]
  levels(bin2.begin) <- tstrsplit(as.character(levels(bin2.begin)), "[][,)]")[[2]]
  levels(bin2.end) <- tstrsplit(as.character(levels(bin2.end)), "[][,)]")[[3]]
  mat[,begin1:=as.integer(as.character(bin1.begin))]
  mat[,end1:=as.integer(as.character(bin1.end))]
  mat[,begin2:=as.integer(as.character(bin2.begin))]
  mat[,end2:=as.integer(as.character(bin2.end))]
  return(mat)
}

#' Build (grouped) signal matrix using normalization data if available
#' @keywords internal
get_signal_matrix = function(cs, resolution=cs@settings$base.res, groups=cs@experiments[,.(name,groupname=name)]) {
  #build signal matrix
  if (resolution >= cs@biases[,max(pos)-min(pos)]) {
    sbins=c(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1])
    signal.bins=unique(cut(c(sbins,head(sbins,n=-1)+resolution/2), sbins,
                           ordered_result=T, right=F, include.lowest=T,dig.lab=12))
    mat=CJ(name=groups[,groupname],bin1=signal.bins,bin2=signal.bins,sorted=F,unique=F)[bin2>=bin1]
    mat[,phi:=0]
  } else {
    #report phi values for each group
    if (resolution!=cs@settings$base.res) {
      refmat=groups[cs@par$signal[,.(name,bin1,bin2,phi)]][
        !is.na(groupname),.(name=groupname,refbin1=bin1,refbin2=bin2,phi)][
          ,.(phi=mean(phi)),by=c("name","refbin1","refbin2")]
      #merge signal to new binning
      sbins=seq(cs@biases[,min(pos)-1],cs@biases[,max(pos)+1+resolution],resolution)
      pos=head(sbins,n=-1)+resolution/2
      bins=unique(data.table(refbin=cut(pos, cs@settings$sbins,
                                        ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                             bin=cut(pos, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)))
      stopifnot(bins[,.N]==length(sbins)-1)
      mat=CJ(name=groups[,groupname],bin1=bins[,unique(bin)],bin2=bins[,unique(bin)],
             sorted=F,unique=F)[bin2>=bin1]
      mat=merge(mat,bins,by.x="bin1",by.y="bin",all.y=T)
      mat=merge(mat,bins,by.x="bin2",by.y="bin",all.y=T, suffixes=c("1","2"))
      mat=merge(mat,refmat,by=c("name","refbin1","refbin2"),all.x=T)
      mat[is.na(phi),phi:=0]
      mat=mat[,.(name,bin1,bin2,phi)]
    } else {
      mat=groups[cs@par$signal[,.(name,bin1,bin2,phi)]][!is.na(groupname),.(name=groupname,bin1,bin2,phi)]
    }
    mat=mat[,.(phi=mean(phi)),keyby=c("name","bin1","bin2")]
  }
  #ggplot(cs@par$signal)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
  #ggplot(mat)+geom_raster(aes(bin1,bin2,fill=phi))+scale_fill_gradient2()
  return(mat)
}


