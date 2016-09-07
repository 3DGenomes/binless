#' @include csnorm.R
NULL

#' Bin a counts data.table into a matrix of a given resolution
#'
#' @param counts data.table as returned by \code{\link{prepare_for_sparse_cs_norm}}
#' @param resolution positive integer.
#' @param b1,b2,e1,e2 Begins and ends of the portion of the data to bin. If NULL replace by min/max value. 
#'
#' @return a data.table representing the binned data
#' @export
#'
#' @examples
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

#' Apply the ICE algorithm to a binned matrix
#'
#' @param csb a CSbinned object
#' @param niterations positive integer. Number of iterations to perform
#'
#' @return a CSbinned object containing the ICEd matrix
#' @export
#'
#' @examples
iterative_normalization = function(raw, niterations=100, namecol="name") {
  binned = foreach (n=raw[,unique(get(namecol))], .combine="rbind") %do% {
    binned = raw[get(namecol)==n&bin1<bin2,.(bin1,bin2,N=observed)]
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
    binned[,c("b1","b2","N"):=list(NULL,NULL,NULL)]
    setnames(binned,"N.weighted",paste0("ice.",niterations))
    binned=binned[bin1<bin2]
    binned[,c(namecol):=n]
    setkeyv(binned,c(namecol,"bin1","bin2"))
  }
}

#' Bin observed and expected counts at a given resolution
#' @keywords internal
#' @export
#' 
csnorm_predict_binned = function(cs, resolution, ncores=1) {
  Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
  Dsets=cs@experiments[,.N]
  stopifnot(Dsets*Kdiag==length(cs@par$beta_diag_centered))
  npoints=100*Kdiag #evaluate spline with 100 equidistant points per basis function
  #bin existing counts and biases
  biases=cs@biases
  counts=cs@counts
  bins=seq(biases[,min(pos)-1],biases[,max(pos)+1+resolution],resolution)
  counts[,c("bin1","bin2"):=list(cut(pos1, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=5),
                                    cut(pos2, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=5))]
  counts[,c("ibin1","ibin2"):=list(as.integer(bin1)-1,as.integer(bin2)-1)]
  biases[,bin:=cut(pos, bins, ordered_result=T, right=F, include.lowest=T,dig.lab=5)]
  biases[,ibin:=as.integer(bin)-1]
  biases[,c("log_nu","log_delta"):=list(cs@par$log_nu,cs@par$log_delta)]
  #split computation across cores
  stepsize=max(2,ceiling(Dsets*length(bins)/(5*ncores)))
  counts[,c("chunk1","chunk2"):=list(ibin1 %/% stepsize, ibin2 %/% stepsize)]
  biases[,chunk:=ibin %/% stepsize]
  chunks=CJ(biases[,min(chunk):max(chunk)],biases[,min(chunk):max(chunk)],cs@experiments[,1:.N])[V1<=V2]
  #run it
  registerDoParallel(cores=ncores)
  binned = foreach (i=chunks[,V1], j=chunks[,V2], k=chunks[,V3], .packages=c("data.table","rstan"), .combine="rbind") %dopar% {
    n=cs@experiments[k,name]
    counts=copy(counts[name==n&chunk1==i&chunk2==j])
    biases1=copy(biases[name==n&chunk==i])
    biases2=copy(biases[name==n&chunk==j])
    if (biases1[,.N]==0 | biases2[,.N]==0) {
      data.table()
    } else {
      coffset=c(biases1[,min(id)-1],biases2[,min(id)-1])
      cidx=t(data.matrix(counts[,.(id1-coffset[1],id2-coffset[2])]))
      boffset=c(biases1[,min(ibin)-1],biases2[,min(ibin)-1])
      data = list( S1=biases1[,.N], S2=biases2[,.N], cutsites1=as.array(biases1[,pos]), cutsites2=as.array(biases2[,pos]),
                   #
                   N=counts[,.N], cidx=cidx,
                   counts_close=as.array(counts[,contact.close]), counts_far=as.array(counts[,contact.far]),
                   counts_up=as.array(counts[,contact.up]), counts_down=as.array(counts[,contact.down]),
                   #
                   B1=biases1[,max(ibin)-min(ibin)+1], B2=biases2[,max(ibin)-min(ibin)+1],
                   bbins1=as.array(biases1[,ibin-boffset[1]]), bbins2=as.array(biases2[,ibin-boffset[2]]),
                   cbins1=as.array(counts[,ibin1-boffset[1]]), cbins2=as.array(counts[,ibin2-boffset[2]]),
                   #
                   Kdiag=Kdiag, npoints=npoints, circularize=cs@settings$circularize,
                   dmin=cs@settings$dmin, dmax=cs@settings$dmax,
                   #
                   eC=cs@par$eC[k],
                   log_nu1=as.array(biases1[,log_nu]), log_nu2=as.array(biases2[,log_nu]),
                   log_delta1=as.array(biases1[,log_delta]), log_delta2=as.array(biases2[,log_delta]),
                   beta_diag_centered=cs@par$beta_diag_centered[((k-1)*Kdiag+1):(k*Kdiag)]
      )
      out <- capture.output(binned <- optimizing(csnorm:::stanmodels$predict_binned,
                                                 data=data, as_vector=F, hessian=F, iter=1, verbose=T, init=0))
      #format output
      so = data.table(melt(binned$par$observed))
      setnames(so,"value","observed")
      so[,expected:=melt(binned$par$expected)$value]
      so[,ncounts:=melt(binned$par$ncounts)$value]
      so[,decaymat:=melt(binned$par$decaymat)$value]
      so=so[ncounts>0]
      so[,c("Var1","Var2"):=list(Var1+boffset[1],Var2+boffset[2])]
      setkey(so,Var1)
      so=biases1[,.(bin1=bin[1]),keyby=ibin][so]
      setkey(so,Var2)
      so=biases2[,.(bin2=bin[1]),keyby=ibin][so]
      so[,c("ibin","i.ibin"):=list(NULL,NULL)]
      so[,name:=n]
      so
    }
  }
  setkey(binned,name,bin1,bin2)
  #remove created entries
  counts[,c("bin1","bin2","ibin1","ibin2","chunk1","chunk2"):=list(NULL,NULL,NULL,NULL,NULL,NULL)]
  biases[,c("bin","ibin","chunk","log_nu","log_delta"):=list(NULL,NULL,NULL,NULL,NULL)]
  return(binned)
}

#' Predict dispersions for a binned matrix
#' @keywords internal
#' @export
#' 
get_dispersions = function(binned, iter=10000) {
  data=list(B=binned[,.N],observed=binned[,observed],expected=binned[,expected],ncounts=binned[,ncounts])
  out <- capture.output(op <- optimizing(stanmodels$dispersions,
                                         data=data, as_vector=F, hessian=F, iter=iter, verbose=T, init=0, init_alpha=1e-5)$par)
  return(op)
}


#' Predict dispersions for a binned matrix
#' @keywords internal
#' @export
#' 
get_dispersions2 = function(binned, iter=10000) {
  data=list(B=binned[,.N],observed=binned[,observed],expected=binned[,expected])
  out <- capture.output(op <- optimizing(stanmodels$dispersions2,
                                         data=data, as_vector=F, hessian=F, iter=iter, verbose=T, init=0, init_alpha=1e-5)$par)
  return(op)
}


#' compute \eqn{p(\Gamma_2>\Gamma_1) = \int_{0}^{+\infty} dx p_{\Gamma_2}(x) \int_{0}^{x} dy p_{\Gamma_1}(y)}
#' @keywords internal
#' 
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

#' compute \eqn{p(\mathcal{N}_2>\mathcal{N}_1) = \int_{-infty}^{+infty} dx p_{\mathcal{N}_2}(x) \int_{-infty}^{x} dy p_{\mathcal{N}_1}(y)}
#' @keywords internal
#' 
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

#' Interaction calling
#' 
#' @param mat
#' @param name character. Name of reference experiment. Used for naming
#'   probability column "prob.gt.name".
#' @param threshold significance threshold, between 0 and 1
#' @param ncores number of cores used for parallelization
#' @param normal.approx integer. Use normal approximation if dispersion is 
#'   larger than this (should be larger than 10)
#'   
#' @return the matrix with called interactions
#' @keywords internal
#' @export
#' 
#' @examples
call_interactions = function(mat, name, threshold=0.95, ncores=1, normal.approx=100) {
  colname=paste0("prob.gt.",name)
  mat[alpha1<normal.approx,c(colname,"detection.type"):=list(
    compute_gamma_overlap(alpha1,beta1,alpha2,beta2,ncores=ncores),"gamma")]
  mat[,c("mean1","sd1"):=list(alpha1/beta1, sqrt(alpha1)/beta1)]
  mat[,c("mean2","sd2"):=list(alpha2/beta2, sqrt(alpha2)/beta2)]
  mat[alpha1>=normal.approx,c(colname,"detection.type"):=list(
    compute_normal_overlap(mean1,sd1, mean2, sd2, ncores=ncores),"normal")]
  mat[,c(colname):=as.numeric(get(colname))]
  mat[,is.significant:=get(colname)>threshold | 1-get(colname)>threshold]
  if (mat[is.na(get(colname)),.N]>0) {
    cat("Warning:",mat[is.na(get(colname)),.N],"out of ",mat[,.N], "calls could not be computed!\n")
  }
  return(mat)
}

#' estimates the values of the count or dispersion required to cross a given threshold
#' 
#' For illustration purposes. Has a border effect at high counts that can be safely ignored.
#'
#' @param observed,expected,dispersion,threshold floats used to generate the plots 
#' @param counts.range scan counts in that range
#' @param disp.range scan dispersion in that range
#' @param compute whether to compute the values that meet the threshold or not. This step can fail if the ranges are set badly.
#'
#' @return a list of three plots
#' @export
#'
#' @examples
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


#' Generate nu and delta genomic biases on evenly spaced points along the genome
#'
#' @param biases data.table.
#' @param beta_nu,beta_delta vectors. spline parameters
#' @param bf_per_kb number of basis functions per kb
#' @param points_per_kb number of evaluation points per kb
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
generate_genomic_biases = function(biases, beta_nu, beta_delta, bf_per_kb=1, points_per_kb=100) {
  begin=biases[,min(pos)]
  end=biases[,max(pos)]
  genome_sz=end-begin
  Krow=round(bf_per_kb*genome_sz/1000)
  S=round(points_per_kb*genome_sz/1000)
  op=optimizing(stanmodels$gen_genomic_biases, data=list(Krow=Krow, S=S, begin=begin, end=end,
                                                         beta_nu=beta_nu, beta_delta=beta_delta),
                as_vector=F, hessian=F, iter=1, verbose=F, init=0)
  dt=as.data.table(op$par)
  setkey(dt, pos)
  return(dt)
}

#' Bin normalized datasets and compute dispersions at that resolution
#' 
#' @param cs CSnorm object, optimized.
#' @param resolution integer. The desired resolution of the matrix.
#' @param ncores integer. The number of cores to parallelize on.
#' @param ice integer. If positive, perform the optional Iterative Correction 
#'   algorithm, useful for comparison purposes. The value determines the number
#'   of iterations.
#' @param dispersion.type integer. Type 1: original dispersion. Type 2: Fit same
#'   dispersion to all bins. Type 3: Fit one dispersion and multiply by ncounts
#'   to get that of a bin.
#' @param verbose
#'   
#' @return A CSnorm object containing an additional CSbinned object in cs@binned
#' @export
#' 
#' @examples
bin_all_datasets = function(cs, resolution=10000, ncores=1, ice=-1, dispersion.type=c(1,2,3), verbose=T) {
  stopifnot(length(dispersion.type)==1 && dispersion.type>=1 && dispersion.type<=3)
  if (get_cs_binned_idx(cs,resolution,dispersion.type,raise=F)>0) {
    stop("Refusing to overwrite already existing matrices at ", resolution/1000, "kb with dispersion type ",dispersion.type,
         ". Use them or remove them from the cs@binned list")
  }
  if (verbose==T) cat("*** build binned matrices for each experiment\n")
  mat=csnorm_predict_binned(cs, resolution, ncores=ncores)
  setkey(mat,name,bin1,bin2)
  if (ice>0) {
    if (verbose==T) cat("*** iterative normalization with ",ice," iterations\n")
    raw=mat[,.(name,bin1,bin2,observed)]
    setkey(raw,name,bin1,bin2)
    iced=iterative_normalization(raw, niterations=ice, namecol="name")
    setkey(iced,name,bin1,bin2)
    mat=merge(mat,iced,all.x=T,all.y=F)
  }
  #write begins/ends
  if (verbose==T) cat("*** write begin/end positions\n")
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
  #store dispersions
  if (verbose==T) cat("*** Dispersion type ",dispersion.type,"\n")
  mat[,dispersion:=switch(dispersion.type,
                          cs@par$alpha,
                          get_dispersions2(mat)$dispersion,
                          get_dispersions(mat)$dispersion)]
  #compute normalized matrices
  mat[,expected.sd:=sqrt(expected+expected^2/dispersion)] #negative binomial SD
  mat[,normalized:=(dispersion+observed)/(dispersion+expected)] #posterior predictive mean and SD
  mat[,normalized.sd:=sqrt(dispersion+observed)/(dispersion+expected)]
  mat[,c("icelike","icelike.sd"):=list(normalized*decaymat,normalized.sd*decaymat)]
  #create CSmatrix and CSbinned object
  csm=new("CSmatrix", mat=mat, group="all", ice=(ice>0), ice.iterations=ifelse(ice>0,ice,NA),
          names=as.character(mat[,unique(name)]), dispersion.fun=deparse(NA))
  csb=new("CSbinned", resolution=resolution, dispersion.type=dispersion.type, grouped=list(csm),
          individual=copy(mat[,.(name,bin1,begin1,end1,bin2,begin2,end2,observed,expected,ncounts,decaymat,dispersion)]))
  cs@binned=append(cs@binned,csb)
  cs
}

#' Group binned matrices of datasets
#'
#' @param cs CSnorm object
#' @param resolution,dispersion.type see \code{\link{bin_all_datasets}}, used to identify the input matrices.
#' @param group The type of grouping to be performed. Any combination of the given arguments is possible.
#' @param dispersion.fun function. How to treat dispersions. Any many-to-one mapping such as sum, min or mean.
#' @inheritParams bin_all_datasets
#'
#' @return CSnorm object
#' @export
#'
#' @examples
group_datasets = function(cs, resolution, dispersion.type,
                          group=c("condition","replicate","enzyme","experiment"), dispersion.fun=sum, ice=-1, verbose=T) {
  #fetch and check inputs
  experiments=cs@experiments
  csbi=get_cs_binned_idx(cs, resolution=resolution, dispersion.type=dispersion.type, raise=T)
  csb=cs@binned[[csbi]]
  group=match.arg(group, several.ok=T)
  dispersion.fun=match.fun(dispersion.fun)
  if (get_cs_matrix_idx(csb, group, dispersion.fun, raise=F)>0) {
    stop("Refusing to overwrite already existing ", group, " group matrices with dispersion function ",deparse(dispersion.fun),
         ". Use them or remove them from the cs@binned[[",csbi,"]]@grouped list and @metadata table")
  }
  #
  if (verbose==T) cat("*** creating groups\n")
  groups=experiments[,.(name,groupno=.GRP,groupname=do.call(paste,mget(group))),by=group][,.(name,groupno,groupname)]
  if (verbose==T) cat("*** merging matrices\n")
  mat=merge(csb@individual,groups,by="name",all=T)[,.(observed=sum(observed),expected=sum(expected),decaymat=mean(decaymat),
                                                      dispersion=dispersion.fun(dispersion)),
                                                   by=c("groupno","groupname","bin1","bin2","begin1","end1","begin2","end2")]
  setnames(mat,"groupname","name")
  setkey(mat,name,bin1,bin2)
  stopifnot(mat[,uniqueN(groupno)==uniqueN(name)]) #otherwise need to change subsequent steps
  mat[,groupno:=NULL]
  #ice
  if (ice>0) {
    cat("*** iterative normalization with ",ice," iterations\n")
    raw=mat[,.(name,bin1,bin2,observed)]
    setkey(raw,name,bin1,bin2)
    iced=iterative_normalization(raw, niterations=ice, namecol="name")
    setkey(iced,name,bin1,bin2)
    mat=merge(mat,iced,all.x=T,all.y=F)
  }
  #compute normalized matrices
  mat[,expected.sd:=sqrt(expected+expected^2/dispersion)] #negative binomial SD
  mat[,normalized:=(dispersion+observed)/(dispersion+expected)] #posterior predictive mean and SD
  mat[,normalized.sd:=sqrt(dispersion+observed)/(dispersion+expected)]
  mat[,c("icelike","icelike.sd"):=list(normalized*decaymat,normalized.sd*decaymat)]
  #store matrices
  csm=new("CSmatrix", mat=mat, group=group, ice=(ice>0), ice.iterations=ifelse(ice>0,ice,NA),
          names=as.character(mat[,unique(name)]), dispersion.fun=deparse(dispersion.fun,nlines=1))
  csb@grouped=append(csb@grouped,list(csm))
  cs@binned[[csbi]]=csb
  return(cs)
}

#' Detect significant interactions wrt expected
#' 
#' @param cs CSnorm object
#' @param resolution,group,dispersion.type,dispersion.fun see
#'   \code{\link{bin_all_datasets}} and \code{\link{group_datasets}}, used to
#'   identify the input matrices.
#' @param threshold significance threshold, between 0 and 1
#' @param ncores number of cores used for parallelization
#' @inheritParams call_interactions
#'   
#' @return the binned matrix with additional information relating to these 
#'   significant interactions
#' @export
#' 
#' @examples
detect_interactions = function(cs, resolution, group, dispersion.type, dispersion.fun, threshold=0.95, ncores=1, normal.approx=100){
  #get CSmat object
  idx1=get_cs_binned_idx(cs, resolution, dispersion.type, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, dispersion.fun, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csm, type="interactions", threshold=threshold, normal.approx=normal.approx, ref="expected", raise=F)>0) {
    stop("Refusing to overwrite this already detected interaction")
  }
  mat=copy(csm@mat[,.(name,bin1,bin2,observed,expected,dispersion)])
  #report gamma parameters
  mat[,c("alpha1","beta1"):=list(dispersion,dispersion)] #posterior
  mat[,c("alpha2","beta2"):=list(alpha1+observed,beta1+expected)] #posterior predictive
  mat=call_interactions(mat,"expected",threshold=threshold,ncores=ncores,normal.approx=normal.approx)
  mat[,c("alpha1","alpha2","beta1","beta2","mean1","mean2","sd1","sd2"):=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)]
  #store back
  csi=new("CSinter", mat=mat, type="interactions", threshold=threshold, normal.approx=normal.approx, ref="expected")
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}


#' Detect significant differences with a reference
#' 
#' @param binned as returned by \code{\link{csnorm_predict_binned}}
#' @param ref character. Detect differences between this group and all others
#' @inheritParams call_interactions
#'   
#' @return the binned matrix with additional information relating to these
#'   significant interactions
#' @export
#' 
#' @examples
detect_differences = function(cs, resolution, group, dispersion.type, dispersion.fun, ref, threshold=0.95, ncores=1, normal.approx=100){
  idx1=get_cs_binned_idx(cs, resolution, dispersion.type, raise=T)
  csb=cs@binned[[idx1]]
  idx2=get_cs_matrix_idx(csb, group, dispersion.fun, raise=T)
  csm=csb@grouped[[idx2]]
  #check if interaction wasn't calculated already
  if (get_cs_interaction_idx(csm, type="differences", threshold=threshold, normal.approx=normal.approx, ref="expected", raise=F)>0) {
    stop("Refusing to overwrite this already detected interaction")
  }
  mat=copy(csm@mat[,.(name,bin1,bin2,observed,expected,dispersion)])
  names=mat[,unique(name)]
  if (!(ref %in% names)) stop("Group name not found! Valid names are: ",deparse(as.character(names)))
  names=names[names!=ref]
  refmat=mat[name==ref,.(bin1,bin2,alpha1=dispersion+observed,beta1=dispersion+expected)]
  setkey(refmat,bin1,bin2)
  qmat = foreach (n=names,.combine="rbind") %do% {
    #merge parameters
    qmat=mat[name==n,.(bin1,bin2,alpha2=dispersion+observed,beta2=dispersion+expected)]
    setkey(qmat,bin1,bin2)
    qmat=merge(refmat,qmat,all=T)
    qmat=qmat[!(is.na(alpha1)|is.na(beta1)|is.na(alpha2)|is.na(beta2))]
    #call interactions
    qmat=call_interactions(qmat,ref,threshold=threshold,ncores=ncores,normal.approx=normal.approx)
    #compute ratio matrix
    qmat[,ratio:=ifelse(alpha2>1,(beta2/beta1)*(alpha1/(alpha2-1)),NA)]
    qmat[,ratio.sd:=ifelse(alpha2>2,(beta2/beta1)*sqrt((alpha1*(alpha1+alpha2-1))/((alpha2-2)*(alpha2-1)^2)),NA)]
    #remove extra columns and add name
    qmat[,c("name","alpha1","alpha2","beta1","beta2","mean1","mean2","sd1","sd2"):=list(n,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)]
    qmat
  }
  #merge everything back, removing extra columns first
  matnames=names(mat)
  if ("is.significant" %in% matnames) {
    mat[,c(matnames[grep(".gt.",matnames)],"detection.type","is.significant"):=list(NULL,NULL,NULL)]
    if ("ratio" %in% matnames) mat[,c("ratio","ratio.sd"):=list(NULL,NULL)]
  }
  mat=merge(mat,qmat,all=T,by=c("name","bin1","bin2"))
  #add interaction to cs object
  csi=new("CSinter", mat=mat, type="differences", threshold=threshold, normal.approx=normal.approx, ref=ref)
  csm@interactions=append(csm@interactions,list(csi))
  csb@grouped[[idx2]]=csm
  cs@binned[[idx1]]=csb
  return(cs)
}
  
