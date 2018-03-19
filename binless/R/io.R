#' @include binless.R
NULL

#' Save normalized CS object, discarding everything but
#' input counts and final normalization parameters
#' 
#' @param cs a normalized CSnorm object
#' @param fname where to save it to
#'   
#' @export
#' 
#' @examples
save_stripped = function(cs, fname) {
  cs@diagnostics=list(params=data.table())
  cs@zeros=data.table()
  strip_interaction = function(csi) {
    if (class(csi) %in% c("CSbsig","CSbdiff")) {
      csi@state=list()
      csi@cts=data.table()
    }
    return(csi)
  }
  strip_group = function(csg) {
    csg@cts=data.table()
    csg@interactions=lapply(csg@interactions,strip_interaction)
    return(csg)
  }
  cs@groups=lapply(cs@groups,strip_group)
  save(cs,file=fname)
}

#' Load normalized CS object saved with save_stripped,
#' rebuild zeros data.table
#' 
#' @param fname file name
#' @param ncores number of CPUs to use to reconstruct zeros
#'   
#' @return A cs object on which signal detection can be performed
#'
#' @export
#' 
#' @examples
load_stripped = function(fname, ncores=1) {
  cat("Loading stripped file\n")
  cs=get(load(fname))
  cat("Filling zeros\n")
  cs@zeros = binless:::get_nzeros(cs, cs@settings$sbins, ncores=ncores)
  cat("Computing cts\n")
  groups=cs@groups
  cs@groups=list()
  cs@groups = foreach (csg=groups) %do% {
    cs=group_datasets(cs, group=csg@group, resolution=csg@resolution, ncores=ncores)
    csg.new=tail(cs@groups,1)[[1]]
    csg.new@interactions = csg@interactions #do not populate each csi@cts
    csg.new
  }
  return(cs)
}


#' Diagnostics plots to monitor convergence of normalization (gaussian
#' approximation)
#' 
#' @param cs a normalized CSnorm object
#' @param start,end the starting and ending steps on which to plot the info. Default is full range.
#'   
#' @return Two plots in a list. The first is for the four log-likelihoods, the
#'   second is for the parameters.
#' @export
#' 
#' @examples
plot_diagnostics = function(cs, start=1, end=Inf) {
  if (is.infinite(end) || cs@diagnostics$params[,max(step)]<end) end = cs@diagnostics$params[,max(step)] 
  plot=ggplot(cs@diagnostics$params[step>=start&step<end+1,.(step,leg,value)])+
    geom_line(aes(step,value))+geom_point(aes(step,value))+facet_wrap(~leg, scales = "free_y")+
    theme(legend.position="bottom")+scale_x_continuous(breaks = seq(start,end,2), minor_breaks = seq(start+1,end,2))
  vars=foreach(var=c("eC","eRJ","eDE","alpha","lambda_iota","lambda_rho","lambda_diag", "lambda1", "lambda2", "eCprime"),
               trans=(c(NA,NA,NA,NA,"log10","log10","log10","log10","log10",NA)),.combine=rbind) %do% get_all_values(cs,var,trans)
  plot2=ggplot(vars[step>=start&step<end+1])+geom_line(aes(step,value,group=as.integer(step)))+
    geom_point(aes(step,value,colour=leg))+facet_wrap(~variable, scales = "free_y")+
    scale_x_continuous(breaks = seq(start,end,2), minor_breaks = seq(start+1,end,2))
  return(list(plot=plot,plot2=plot2))
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
  p=ggplot(binned, aes(begin1,begin2, fill=log(observed)))+geom_raster()+
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

#' make plot of binless matrix
#'
#' @param mat the binless matrix
#' @param upper,lower the name of the column to plot in upper (resp. lower)
#'  triangle. Default: phi. Formulae are allowed (ggplot2 aes syntax)
#' @param facet what to facet on. Default: name. Pass a character vector for multiple facets
#' @param limits set to a pair of values to enforce plot and colour scale to be within these values
#'
#' @return
#' @export
#'
#' @examples
plot_binless_matrix = function(mat, upper="binless", lower="observed", trans="log10", facet="name", limits=NULL, label=NA) {
  data=copy(mat)
  if (!("begin1" %in% names(data) && "begin2" %in% names(data))) {
    if (data[,is.factor(bin1)] && data[,is.factor(bin2)]) {
      data = add_bin_begin_and_end(data)
      bin1="begin1"
      bin2="begin2"
      just=1
    } else {
      bin1="bin1"
      bin2="bin2"
      just=0.5
    }
  } else {
    bin1="begin1"
    bin2="begin2"
    just=1
  }
  #plot data
  p=ggplot(data)+geom_raster(aes_string(bin1,bin2,fill=upper),hjust=just,vjust=just)+
    geom_raster(aes_string(bin2,bin1,fill=lower),hjust=just,vjust=just)
  if (all(facet %in% names(data))) p = p + facet_wrap(facet)
  #set colour scale
  p=p+scale_fill_gradient2(low=muted("blue"),high=muted("red"),na.value = "white", limits=limits, oob=squish, trans=trans)
  #set visual appearance
  p=p+coord_fixed()+theme_minimal()+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.spacing=unit(0,"cm"), axis.title=element_blank())+
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  if (!is.na(label)) p = p + labs(fill=label)
  print(p)
  invisible(p)
}

#' Plot a binless signal matrix
#' 
#' This function is to ensure a consistent representation of signal matrices
#' Plots the signal in log2 scale (unit is log2 fold change) and caps the maximum
#' signal to log2FC = 3.
#'
#' @param mat the binless matrix data.table
#'
#' @return
#' @export
#'
#' @examples
plot_binless_signal_matrix = function(mat) {
  plot_binless_matrix(mat,upper="log2(signal)", lower="log2(signal)", trans="identity", limits=c(-3,3), label="log2 FC")
}

#' Plot a binless difference matrix
#' 
#' Plots the difference in log2 scale (unit is log2 fold change) and caps the maximum
#' difference (in absolute value) to log2FC = 3.
#'
#' @param mat the binless matrix data.table
#'
#' @return
#' @export
#'
#' @examples
plot_binless_difference_matrix = function(mat) {
  plot_binless_matrix(mat,upper="log2(difference)", lower="log2(difference)", trans="identity", limits=c(-3,3), label="log2 FC")
}

