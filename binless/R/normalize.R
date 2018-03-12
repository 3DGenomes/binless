#' @include binless.R
NULL

#' count number of zeros in a given cut site, distance bin and signal bin
#' 
#' @return a data table, keyed on name,dbin,id1,pos1,bin1,bin2,dir,cat with the following columns
#' - name: name of the dataset
#' - id1, pos1, bin1: coordinates of the cut site
#' - bin2, dbin: signal/distance bin in which we looked at the intersections
#' - dir: fwd (rev) for contacts with downstream (resp upstream) cut-sites
#' - cat: whether we consider contacts on the left (contact L) or on the right (contact R) of this cut site
#' - ncross: number of cut site intersections (crossings) in this signal/distance/direction/category bin.
#'  We discard anything below cs@settings$dmin. Note that sum(ncross) is four times the total number of crossings (per dataset).
#' - nnz: number of non-zero contacts in this signal/distance/direction bin (max 2 per crossing).
#'  Note that sum(nnz) is twice the number of nonzeros (per dataset)
#' - nzero: number of zeros in this signal/distance/direction bin. We have nzero = 2*ncross - nnz.
#'  Note that sum(nzero) is twice the number of zeros (per dataset)
#'  For speed purposes downstream, we only return the entries where nzero>0, which corresponds to most of the entries anyway.
#'  We therefore have sum(nnz+nzero) equal to approximately twice the number of detectable counts, e.g. 8x the number of crossings
#'  
#' @keywords internal
#' 
get_nzeros = function(cs, sbins, ncores=1) {
  stopifnot(cs@counts[id1>=id2,.N]==0)
  #positive counts: group per cut site and signal / distance bin
  cts=melt(cs@counts[,.(name,pos1,pos2,distance,contact.close,contact.down,contact.far,contact.up)],
           id.vars=c("name","pos1","pos2","distance"))[value>0]
  cts[,c("bin1","bin2","dbin"):=
        list(cut(pos1, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
             cut(pos2, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
             cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12))]
  cts=rbind(cts[,.(name,pos1,bin1,bin2,dbin,variable)][
              ,.(nnz=.N,dir="fwd"),keyby=c("name","pos1","bin1","bin2","dbin","variable")],
            cts[,.(name,pos1=pos2,bin1=bin2,bin2=bin1,dbin,variable)][
              ,.(nnz=.N,dir="rev"),keyby=c("name","pos1","bin1","bin2","dbin","variable")])
  cts[,cat:=ifelse(dir=="fwd",ifelse(variable %in% c("contact.up","contact.far"), "contact L", "contact R"),
                              ifelse(variable %in% c("contact.up","contact.close"), "contact L", "contact R"))]
  cts=cts[,.(nnz=sum(nnz)),keyby=c("name","pos1","bin1","bin2","dbin","cat","dir")]
  #Count the number of crossings per distance bin
  #looping over IDs avoids building NxN matrix
  biases=cs@biases[,.(name,id,pos)]
  biases[,bin:=cut(pos, sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
  chunksize=cs@biases[,ceiling(.N/(10*ncores))]
  nchunks=cs@biases[,ceiling(.N/chunksize)]
  registerDoParallel(cores=ncores)
  crossings = foreach(chunk=1:nchunks, .combine=rbind) %dopar% {
    bs=biases[((chunk-1)*chunksize+1):min(.N,chunk*chunksize)]
    foreach(i=bs[,id], n=bs[,name], p=bs[,pos], b=bs[,bin], .combine=rbind) %do% {
      crossings = biases[name==n&pos!=p,.(name,pos2=pos,bin2=bin,distance=abs(pos-p))]
      if (cs@settings$circularize>0)  crossings[,distance:=pmin(distance,cs@settings$circularize+1-distance)]
      crossings[,c("dbin","dir"):=list(
        cut(distance,cs@settings$dbins,ordered_result=T,right=F,include.lowest=T,dig.lab=12),
        ifelse(pos2>p,"fwd","rev"))]
      crossings[distance>=cs@settings$dmin,.(id1=i, pos1=p,bin1=b,ncross=.N),by=c("name","bin2","dbin","dir")]
    }
  }
  stopImplicitCluster()
  #now merge with positive counts and deduce number of zeros
  zeros = rbind(
    merge(crossings,cts[cat=="contact L"],by=c("name","pos1","bin1","bin2","dbin","dir"),all=T)[
      ,.(name,id1,pos1,bin1,bin2,dbin,dir,cat="contact L",ncross,nnz)],
    merge(crossings,cts[cat=="contact R"],by=c("name","pos1","bin1","bin2","dbin","dir"),all=T)[
      ,.(name,id1,pos1,bin1,bin2,dbin,dir,cat="contact R",ncross,nnz)])
  zeros[is.na(nnz),nnz:=0]
  zeros[,nzero:=2*ncross-nnz]
  stopifnot(zeros[is.na(ncross),.N==0])
  stopifnot(zeros[nzero<0,.N==0])
  zeros=zeros[nzero>0]
  setkey(zeros,name,dbin,id1,pos1,bin1,bin2,dir,cat)
  return(zeros)
}

#' Cleanup a CSnorm object, store settings and populate it with initial guesses of all required parameters
#' @keywords internal
#' 
fresh_start = function(cs, bf_per_kb=50, bf_per_decade=10, bins_per_bf=10, base.res=5000,
                       bg.steps=5, iter=100, fit.signal=T, verbose=T, ncounts=100000, init.dispersion=1,
                       tol=1e-1, ncores=1, fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA) {
    #fresh start
    cs@par=list() #in case we have a weird object
    cs@groups=list()
    cs@diagnostics=list()
    #add settings
    if (fit.signal==F) base.res = cs@biases[,max(pos)-min(pos)]+2
    cs@settings = c(cs@settings[c("circularize","dmin","dmax","qmin","dfuse")],
                    list(bf_per_kb=bf_per_kb, bf_per_decade=bf_per_decade, bins_per_bf=bins_per_bf, base.res=base.res,
                         bg.steps=bg.steps, iter=iter, init.dispersion=init.dispersion, tol=tol,
                         fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                         fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at))
    cs@settings$Kdiag=round((log10(cs@settings$dmax)-log10(cs@settings$dmin))*cs@settings$bf_per_decade)
    cs@settings$Krow=round(cs@biases[,(max(pos)-min(pos))/1000*cs@settings$bf_per_kb])
    stepsz=1/(cs@settings$bins_per_bf*cs@settings$bf_per_decade)
    cs@settings$dbins=10**seq(log10(cs@settings$dmin-1),log10(cs@settings$dmax+1)+stepsz,stepsz)
    #initial guess
    if (verbose==T) cat("No initial guess provided\n")
    decay=CJ(name=cs@experiments[,name],dist=head(cs@settings$dbins,n=length(cs@settings$dbins)-1)*10**(stepsz/2))
    decay[,dbin:=cut(dist, cs@settings$dbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)]
    decay[,c("dist","log_decay"):=list(NULL,0)]
    cs@par=list(eC=array(0,cs@experiments[,.N]), eRJ=array(0,cs@experiments[,.N]), eDE=array(0,cs@experiments[,.N]), alpha=init.dispersion,
                log_iota=array(0,cs@biases[,.N]), log_rho=array(0,cs@biases[,.N]),
                decay=decay, log_decay=0, tol_genomic=.1, tol_decay=.1, tol_disp=.1, tol_signal=1)
    #prepare signal matrix
    if (fit.signal==T) {
      if(verbose==T) cat("Preparing for signal estimation\n")
      stuff = binless:::prepare_first_signal_estimation(cs@biases, cs@experiments[,name], base.res)
      cs@par$signal=stuff$signal
      cs@par$beta.phi=stuff$signal[,beta]
      cs@par$lambda1=array(dim=cs@experiments[,.N])
      cs@settings$sbins=stuff$sbins
    } else {
      cs@settings$sbins=cs@biases[,c(min(pos)-1,max(pos)+1)]
      cs@par$signal=cs@biases[,.(phi=0,beta=0,bin1=cut(pos[1], cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12),
                                 bin2=cut(pos[1], cs@settings$sbins, ordered_result=T, right=F, include.lowest=T,dig.lab=12)),by=name]
      cs@par$beta.phi=cs@par$signal[,beta]
      cs@par$lambda1=NA
      setkey(cs@par$signal,name,bin1,bin2)
    }
    #get number of zeros along cut sites and decay
    if(verbose==T) cat("Counting zeros\n")
    cs@zeros = binless:::get_nzeros(cs, cs@settings$sbins, ncores=ncores)
    #set initial guess for exposures, decay and biases
    if(verbose==T) cat("Initial guess: residuals\n")
    cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
    if(verbose==T) cat("Initial guess: exposures\n")
    cs = binless:::initial_guess_exposures(cs, cts.common)
    if(verbose==T) cat("Initial guess: decay\n")
    cs = binless:::initial_guess_decay(cs, cts.common)
    if(verbose==T) cat("Initial guess: biases\n")
    cs = binless:::initial_guess_genomic(cs, cts.common)
    return(cs)
}

#' update diagnostics data table or create it if not existing
#'
#' @return
#' @keywords internal
#'
#' @examples
update_diagnostics = function(cs, step, leg, runtime) {
  params=data.table(step=step,leg=leg,value=cs@par$value,runtime=runtime)
  tmp=as.data.table(lapply(cs@par,list))
  #remove entries that are too heavy and redundant
  if ("biases" %in% names(tmp)) tmp[,biases:=NULL]
  #merge with previous
  params=cbind(params,tmp)
  if (is.data.table(cs@diagnostics$params)) params=rbind(cs@diagnostics$params,params,fill=T)
  return(params)
}


#' Returns one param's values during optimization
#'
#' @keywords internal
#' @export
#'
#' @examples
get_all_values = function(cs, param, trans) {
  #get value in tmp as vector of lists, remove NULL lists
  legs=c("expo","bias","decay","signal","disp")
  if (!(param %in% names(cs@diagnostics$param))) return(data.table())
  values=cs@diagnostics$params[,.(step,leg=ordered(leg,legs),tmp=get(param))][!sapply(tmp,is.null)]
  values[,step:=step+((unclass(leg)-1)%%10)/10]
  #melt it
  melted=as.data.table(values[,melt(tmp)])
  if ("Var1" %in% names(melted)) {
    if ("Var2" %in% names(melted)) {
      melted[,variable:=paste0(param,".",Var1,".",Var2)]
    } else {
      melted[,variable:=paste0(param,".",Var1)]
    }
  } else {
    melted[,variable:=param]
  }
  #merge it back
  values[,L1:=.I]
  melted=merge(values,melted,by="L1")[,.(variable,step,leg,value)]
  #transform if necessary
  if (!is.na(trans)) {
    if (trans=="log") {
      melted[,c("variable","value"):=list(paste(variable,"(log)"),log(value))]
    } else if (trans=="log10") {
      melted[,c("variable","value"):=list(paste(variable,"(log10)"),log10(value))]
    } else if (trans=="exp") {
      melted[,c("variable","value"):=list(paste(variable,"(exp)"),exp(value))]
    } else {
      stop("unsupported transformation ",trans)
    }
  }
  #return
  setkey(melted,variable,step,leg)
  melted
}

#' Check whether a normalization has converged
#' @export
#' 
has_converged = function(cs) {
  #return FALSE if legs have changed, and require at least 2 steps
  params=cs@diagnostics$params
  laststep=params[,step[.N]]
  if (laststep<=2) return(FALSE)
  if (!setequal(params[step==laststep,.(leg)],params[step==laststep-1,.(leg)])) return(FALSE)
  #check all legs present in the last step
  #check if all quantities directly involved in building the expected matrix have converged
  rel.precision=function(name,fn=identity){merge(params[step==laststep,.(leg,get(name))],
                                                 params[step==laststep-1,.(leg,get(name))],by="leg")[
                                                   leg==leg[.N], max( abs(fn(V2.x[[1]])-fn(V2.y[[1]])) ) / ( max(fn(V2.x[[1]]))-min(fn(V2.x[[1]])) ) ]}
  conv.log_iota = rel.precision("log_iota")
  conv.log_rho = rel.precision("log_rho")
  conv.log_decay = rel.precision("log_decay")
  conv.signal = max(rel.precision("beta.phi"))
  #cat(" conv.log_iota ", conv.log_iota,
  #    " conv.log_rho ", conv.log_rho, " conv.log_decay ", conv.log_decay, " conv.phi ", conv.phi, "\n")
  cat(" relative precision for this iteration: iota ", conv.log_iota,
      " rho ", conv.log_rho, " decay ", conv.log_decay, " signal ", conv.signal, "\n")
  conv.param = all(c(conv.log_iota,conv.log_rho,
                     conv.log_decay,conv.signal)<cs@settings$tol, na.rm=T)
  return(conv.param)
}

#' Check whether a normalization has converged
#'
#' @keywords internal
#' 
get_residuals = function(cts.common, viewpoint) {
  a=cts.common[bin1==viewpoint,.(signal=sum(exp(phi)*nobs),
                               decay=sum(exp(log_decay)*nobs),
                               bias=sum(exp(log_bias)*nobs),
                               mean=sum(exp(lmu.nosig+phi)*nobs),
                               ncounts=sum(nobs),
                               count=sum(count*nobs)),
             keyby=c("name","bin2")]
  setnames(a,"bin2","bin")
  a[bin==viewpoint,c("signal","decay","bias","mean","ncounts","count"):=list(signal/2,decay/2,bias/2,mean/2,ncounts/2,count/2)]
  return(a)
}

#' Binless normalization
#' 
#' @param cs CSnorm object as returned by \code{\link{merge_cs_norm_datasets}}
#' @param init boolean (default is FALSE). Whether to continue an existing
#'  normalization or to start anew (default).
#' @param bf_per_kb positive numeric. Number of cubic spline basis functions per
#'   kilobase, for genomic bias estimates. Small values make the optimization 
#'   easy, but makes the genomic biases stiffer.
#' @param bf_per_decade positive numeric. Number of cubic spline basis functions
#'   per distance decade (in bases), for diagonal decay. Default parameter 
#'   should suffice.
#' @param bins_per_bf positive integer. Number of distance bins to split basis 
#'   functions into. Must be sufficiently small so that the diagonal decay is 
#'   approximately constant in that bin.
#' @param lambdas positive numeric. Length scales to try out as initial
#'   condition.
#' @param ngibbs positive integer. Number of gibbs sampling iterations.
#' @param iter positive integer. Maximum number of optimization steps for each leg.
#' @param fit.signal boolean. Set to FALSE only for diagnostics.
#' @param verbose Display progress if TRUE
#' @param init.dispersion positive numeric. Value of the dispersion to use initially.
#' @param tol positive numeric (default 1e-2). Convergence tolerance on relative changes in the computed biases.
#' @param ncores positive integer (default 1). Number of cores to use.
#' @param fix.lambda1 whether to set lambda1 to a given value, or to estimate it
#' @param fix.lambda1.at if fix.lambda1==T, the approximate value where it is meant to be fixed. Might move a bit because
#'   of the positivity and degeneracy constraints.
#'   
#' @return A csnorm object
#' @export
#' 
#' @examples
#' 
normalize_binless = function(cs, restart=F, bf_per_kb=50, bf_per_decade=10, bins_per_bf=10, base.res=5000,
                     ngibbs = 15, bg.steps=5, iter=100, fit.signal=T,
                     verbose=T, ncounts=100000, init.dispersion=1,
                     tol=1e-1, ncores=1, fix.lambda1=F, fix.lambda1.at=NA, fix.lambda2=F, fix.lambda2.at=NA) {
  #basic checks
  stopifnot( (cs@settings$circularize==-1 && cs@counts[,max(distance)]<=cs@biases[,max(pos)-min(pos)]) |
               (cs@settings$circularize>=0 && cs@counts[,max(distance)]<=cs@settings$circularize/2))
  #
  if (verbose==T) cat("Normalization with fast approximation and performance iteration\n")
  setkey(cs@biases, id, name)
  setkey(cs@counts, id1, id2, name)
  if (restart==F) {
    #fresh start
    cs = fresh_start(cs, bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, base.res = base.res,
                     bg.steps = bg.steps, iter = iter, fit.signal = fit.signal,
                     verbose = verbose, ncounts = ncounts, init.dispersion = init.dispersion,
                     tol = tol, ncores = ncores, fix.lambda1 = fix.lambda1, fix.lambda1.at = fix.lambda1.at,
                     fix.lambda2 = fix.lambda2, fix.lambda2.at = fix.lambda2.at)
    laststep=0
  } else {
    if (verbose==T) cat("Continuing already started normalization with its original settings\n")
    laststep = cs@diagnostics$params[,max(step)]
    cs@groups=list()
  }
  #
  if(verbose==T) cat("Subsampling counts for dispersion\n")
  subcounts = binless:::subsample_counts(cs, ncounts)
  subcounts.weight = merge(cs@zeros[,.(nc=sum(ncross)/4),by=name],subcounts[,.(ns=.N),keyby=name],by="name")[,.(name,wt=nc/ns)]
  #gibbs sampling
  if (ngibbs==0) return(cs)
  for (i in (laststep + 1:ngibbs)) {
    if (verbose==T) cat("\n### Iteration",i,"\n")
    #
    #compute residuals once for this round
    if(verbose==T) cat(" Residuals\n")
    cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
    residuals = binless:::get_residuals(cts.common, viewpoint = cts.common[,min(bin1)])
    residuals[,step:=i]
    cs@diagnostics$residuals = rbind(cs@diagnostics$residuals, residuals)
    #
    #fit exposures
    a=system.time(cs <- binless:::gauss_exposures(cs, cts.common, verbose=verbose))
    cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="expo", runtime=a[1]+a[4])
    #
    #fit iota and rho
    constrain.bias = fit.signal==T && i <= cs@settings$bg.steps+1
    a=system.time(cs <- binless:::gauss_genomic(cs, cts.common, verbose=verbose, constrain=constrain.bias))
    cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="bias", runtime=a[1]+a[4])
    #
    if (fit.signal==T && i > cs@settings$bg.steps) {
      if(verbose==T) cat(" Residuals\n")
      cts.common = binless:::gauss_common_muhat_mean(cs, cs@zeros, cs@settings$sbins)
      #fit signal using sparse fused lasso
      a=system.time(cs <- binless:::gauss_signal(cs, cts.common, verbose=verbose, ncores=ncores,
                                                       fix.lambda1=cs@settings$fix.lambda1,
                                                       fix.lambda1.at=cs@settings$fix.lambda1.at,
                                                       fix.lambda2=cs@settings$fix.lambda2,
                                                       fix.lambda2.at=cs@settings$fix.lambda2.at))
      cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="signal", runtime=a[1]+a[4])
    } else {
      #fit diagonal decay
      a=system.time(cs <- binless:::gauss_decay(cs, cts.common, verbose=verbose))
      cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="decay", runtime=a[1]+a[4])
    }
    #
    #fit dispersion
    a=system.time(cs <- binless:::gauss_dispersion(cs, counts=subcounts, weight=subcounts.weight, verbose=verbose))
    cs@diagnostics$params = binless:::update_diagnostics(cs, step=i, leg="disp", runtime=a[1]+a[4])
    #
    #check for convergence
    if (has_converged(cs)) {
      if (fit.signal == T && i <= cs@settings$bg.steps) {
        cat("Background has converged, fitting signal\n")
        cs@settings$bg.steps = i #compute signal at next step
      } else {
        if (verbose==T) {
          cat("Normalization has converged\n")
        }
        break
      }
    }
  }
  if (verbose==T) cat("Done\n")
  return(cs)
}


