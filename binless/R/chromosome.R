#' @include binless.R
NULL


#' Find optimized parameters and normalize whole chromosome using binless
#'
#' This procedure combines \code{\link{normalize_binless}} and
#' \code{\link{fast_binless}} to propose a chromosome-wide normalization of a
#' dataset, or a difference estimation between two datasets. First, submatrices
#' are generated along the diagonal, and some heuristics (see paper) are
#' computed to select the top 10 (by default) which will likely contain a signal
#' (or a difference). These matrices are subsequently normalized with
#' normalize_binless, followed by signal (difference) detection. The submatrix
#' with the smallest lambda2 is used to provide inputs to fast_binless, which
#' will be computed on the whole dataset.
#'
#' @param sample The path to a CSdata object, stored as RData (see
#'   preprocessing.R), to perform binless signal detection
#' @param sample2 The path to a second CSdata object, in which case binless
#'   differences will be computed
#' @param zoom_region If "full" (default), will normalize complete dataset. If a
#'   vector c(start,end) is passed, will restrict normalization to this portion.
#' @param ncores The number of cores to parallelize parameter estimation
#' @param base.res The base resolution of the normalization
#' @param far.cutoff distance at which we should stop computing signal. Used for
#'   final lambda1 estimation, to lower memory footprint, and to write matrices.
#'   Default is 10e6 (10Mb)
#' @param out_prefix how output files will be named. Can include a path
#' @param write.intermediate Whether to write fast binless output to disk, in
#'   order to save some memory. Only tested in signal calculations
#' @param write.logfile Output of individual optimized_binless calculations will
#'   be written there
#' @param max.submatrices Run optimized_binless on that many submatrices, at
#'   most.
#' @param submat.width.factor Submatrices will have that many bins (multiply by
#'   base.res to get size). Default is 120
#' @param submat.step.factor Submatrices will be extracted every so many bins.
#'   Default is 20
#'
#' @inheritParams normalize_binless
#' @inheritParams fast_binless
#'
#' @return Nothing will be returned. Will write the following files, all
#'   starting with out_prefix. _submatrix_parameters.tsv is a tab-separated file
#'   containing the normalization parameters of all submatrices. _raw_band.RData
#'   is the raw data up to far.cutoff, as a sparseMatrix object;
#'   _binless_band.RData is the binless matrix, up to far.cutoff;
#'   _log10signal.RData is the binless log10(signal) matrix. When differences
#'   are computed, an object is created for each of the reference (first
#'   dataset, _ref in filename) or the target (second dataset, _target in
#'   filename), and two additional matrices are produced: _log10diff.RData
#'   contains the unthresholded difference between target and reference, and
#'   _log10diff_signif.RData uses the threshold provided by the optimization
#'   procedure. It is however recommended to look at both matrices, since this
#'   threshold's estimation is not very robust
#' @export
#'
#' @examples " "
chromosome_binless = function(sample, sample2=NA, zoom_region="full", ncores=16, base.res=5000, far.cutoff=10e6,
                              out_prefix="chromosome_binless", write.intermediate=F, write.logfile=F, max.submatrices=10,
                              submat.width.factor=120, submat.step.factor=20,
                              #common args
                              bg.steps=5, tol_val = 2e-1,
                              #normalize_binless args
                              bf_per_kb=50, bf_per_decade=10, bins_per_bf=10,
                              ngibbs = 25, iter=100, init.dispersion=.1, nrows.dispersion = 100,
                              min.lambda2=.1, fix.lambda1=F, fix.lambda1.at=NA,
                              fix.lambda2=T, fix.lambda2.at=2.5,
                              #fast_binless args
                              nouter = 25, free_decay = 10000) {
  
  ### deduced parameters
  calc.differences = !is.na(sample2)
  maxdiag = far.cutoff/base.res
  submat_step = submat.step.factor*base.res
  submat_step_bin=floor(submat_step/base.res)
  submat_length = submat.width.factor*base.res
  
  ### put csd file(s) in cs object
  start.time = Sys.time()
  if (calc.differences==T) {
      allcsd = list(get(load(sample)),get(load(sample2)))
  } else {
      allcsd = list(get(load(sample)))
  }
  cs=merge_cs_norm_datasets(allcsd, different.decays="none")
  
  if (!(length(zoom_region)==1 && zoom_region == "full")) {
    cs = zoom_csnorm(cs, zoom_region[1], zoom_region[2])
  }
  
  end.time = Sys.time()
  time.prepare =  end.time - start.time
  cat("*** Prepare data time=",time.prepare,"\n")
  
  
  ### create zooms for each submatrix
  start.time = Sys.time()
  
  start_positions = seq.int(from=cs@biases[,min(pos)],to=cs@biases[,max(pos)]+1-submat_length, by=submat_step)
  cs_arr = foreach(start_pos=start_positions) %do% zoom_csnorm(cs,start_pos,start_pos+submat_length)
  
  ### get standard deviation of counts along each submatrix of size submat_length
  ### keep submatrices which have more than average counts
  if (write.logfile==T) {
      cl<-makeCluster(ncores, outfile=paste0(out_prefix,"_cluster.log"))
  } else {
      cl<-makeCluster(ncores)
  }
  registerDoParallel(cl)
  counts_stats=foreach (csa=cs_arr, st=start_positions, .errorhandling = 'remove',
              .packages=c("foreach","data.table","binless"), .combine=rbind) %dopar% {
      matsub = binless:::bin_data(csa,resolution=base.res)
      matsub[,rescaled:=observed/sum(observed),by=name]
      st=matsub[,.(start=st,counts=sum(observed),std=sd(rescaled)),by=name]
      st
  }
  counts_stats[,counts:=counts/sum(counts),by=name]
  counts_stats=counts_stats[,.(counts=mean(counts),std=mean(std)),by=start]
  counts_stats[,is.valid:=counts>mean(counts)]
  cs_arr = cs_arr[counts_stats[,is.valid==T]]
  start_positions = start_positions[counts_stats[,is.valid==T]]
  counts_stats = counts_stats[is.valid==T,.(start,counts,std)]
  
  ### compute directionality index (for signal) or sum of differences on each submatrix 
  if (calc.differences == F ) {
    #signal: compute directionality index on all submatrices
    di=foreach (i=seq_len(counts_stats[,.N]), csa=cs_arr, .errorhandling = 'remove',
                .packages=c("foreach","data.table","binless"), .combine=rbind) %dopar% {
      matsub = binless:::bin_data(csa,resolution=base.res)
      sddi=mean(matsub[,di_index(observed,floor(submat_length/base.res/10),floor((.N*2)**0.5)),by=name]$V1)
      st=cbind(counts_stats[i],data.table(std.directionality=sddi))
      st
    }
    di=di[order(-std.directionality)]
  } else {
      #differences: compute sum of absolute differences larger than 0.1 on all submatrices
    di=foreach (i=seq_len(counts_stats[,.N]), csa=cs_arr, .errorhandling = 'remove',
                .packages=c("foreach","data.table","binless"), .combine=rbind) %dopar% {
       matsub = binless:::bin_data(csa,resolution=base.res)
       alpha=0.1
       lam2=2.5
       lam1=0
       out=binless:::fast_binless(matsub, matsub[,nlevels(bin1)],
                                  alpha, lam2, lam1, nouter, tol_val,
                                  bg.steps, free_decay, compute_patchnos = F,
                                  csv_out="", maxdiag=maxdiag)
       ret=binless:::fast_binless_difference(out,1,alpha,lam2,lam1,tol_val)
       ret=as.data.table(ret)
       st=cbind(counts_stats[i],data.table(diffsum=ret[name==name[.N]&abs(difference)>0.1,sum(abs(difference))]))
       st
    }
    di=di[order(-diffsum)]
  }
  
  ### Keep at most 10 non-overlapping submatrices
  di[,bin:=cut(start,seq.int(min(start),max(start)+submat_length+1,by=submat_length),include.lowest = T, ordered_result = T)]
  di[,is.selected:=start==start[1],by=bin]
  di[,is.selected:=is.selected&cumsum(is.selected)<=max.submatrices]
  cs_arr = cs_arr[di[,is.selected==T]]
  start_positions = start_positions[di[,is.selected==T]]
  #save(start_positions, cs_arr, file="inputs.RData")
  
  ### run parallel normalizations and signal (difference) detections
  ### and extract sensible values for lambda1, lambda2 and alpha
  all_chunks_norm = foreach (start_pos=start_positions, csb=cs_arr, #.errorhandling = 'remove',
              .packages=c("foreach","data.table","binless"), .combine=rbind) %dopar% {
        matsub = binless:::bin_data(csb,resolution=base.res)
        #csb = cs_arr[[c]]
        csb <- normalize_binless(csb, ngibbs = ngibbs, ncores = 2, base.res = base.res, bg.steps = bg.steps, tol = tol_val,
                                 bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter,
                                 min.lambda2 = min.lambda2, init.dispersion=init.dispersion, nrows.dispersion = nrows.dispersion,
                                 fix.lambda1=fix.lambda1, fix.lambda1.at=fix.lambda1.at,
                                 fix.lambda2=fix.lambda2, fix.lambda2.at=fix.lambda2.at)
        
        #save(csb,file=paste0(format(as.integer(start_pos),scientific = FALSE),"-",format(as.integer(start_pos+submat_length),scientific = FALSE),'.RData'))
        if (has_converged(csb, "signal")) {
          csb=bin_all_datasets(csb, ncores = 2)
          if (calc.differences==F) {
              csb=detect_binless_interactions(csb, ncores = 2)
          } else {
              csb=detect_binless_differences(csb, cs@experiments[1,name], ncores = 2)
          }
          warnings()
          #save(csb,file=paste0(format(as.integer(start_pos),scientific = FALSE),"-",format(as.integer(start_pos+submat_length),scientific = FALSE),'.RData'))
          dt=data.table(start=start_pos,end=start_pos+submat_length,
                     alpha=csb@par$alpha[1],
                     lambda2=csb@groups[[1]]@interactions[[1]]@par$lambda2[[1]],
                     lambda1=csb@groups[[1]]@interactions[[1]]@par$lambda1[[1]])
        } else {
          dt=data.table(start=start_pos,end=start_pos+submat_length,
                        alpha=NA,lambda2=NA,lambda1=NA)
        }
        dt
  }
  stopCluster(cl)
  #saveRDS(all_chunks_norm,file=paste0(out_prefix,"_submatrix_parameters.rda"))
  fwrite(all_chunks_norm,file=paste0(out_prefix,"_submatrix_parameters.tsv"),sep="\t")
  prmatrix(all_chunks_norm)
  all_chunks_norm = all_chunks_norm[!is.na(alpha),]
  all_chunks_norm = all_chunks_norm[alpha<100,] #Filter crazy alpha
  all_chunks_norm = all_chunks_norm[!is.na(lambda1),]
  all_chunks_norm = all_chunks_norm[!is.na(lambda2),]
  if(NCOL(all_chunks_norm)>1) {
    minlamidx=as.numeric(which(all_chunks_norm == min(all_chunks_norm$lambda2), arr.ind = TRUE)[1])
    lam1.final=all_chunks_norm[minlamidx]$lambda1
    lam2=all_chunks_norm[minlamidx]$lambda2
    alpha=all_chunks_norm[minlamidx]$alpha
  } else {
    lam1.final=all_chunks_norm$lambda1
    lam2=all_chunks_norm$lambda2
    alpha=all_chunks_norm$alpha
  }
  end.time = Sys.time()
  time.norm = diff(c(start.time, end.time))
  units(time.norm) = "mins"
  cat("*** Optimized binless time=",time.norm,"mins\n")
  
  
  ### Run fast binless on whole chromosome, first binning the whole chromosome
  start.time = Sys.time()
  lam1=0 # do a posteriori thresholding
  cat("*** Obtained lambda1=",lam1.final,"\n")
  cat("*** Using lambda1=",lam1,"\n")
  cat("*** Using lambda2=",lam2,"\n")
  cat("*** Using alpha=",alpha,"\n")
  mat=binless:::bin_data(cs,resolution=base.res)
  
  
  cat("fast binless\n")
  
  if (write.intermediate==T) {
      csvdata=paste0(out_prefix,"_tmp.csv.gz")
      out=binless:::fast_binless(mat, mat[,nlevels(bin1)], alpha, lam2, lam1,
                                 nouter, tol_val, bg.steps, free_decay, compute_patchnos = F,
                                 csv_out=csvdata, maxdiag=maxdiag)
      df=read.csv(csvdata)
      mat = mat[,c('name','bin1','pos1','bin2','pos2','distance','observed','nobs')]
      mat = cbind(mat[,c('name','bin1','pos1','bin2','pos2','distance','observed','nobs')],df)
      out$mat=mat
      rm(df)
      unlink(csvdata)
  } else {
      out=binless:::fast_binless(mat, mat[,nlevels(bin1)], alpha, lam2, lam1,
                                 nouter, tol_val, bg.steps, free_decay, compute_patchnos = F,
                                 csv_out="", maxdiag=maxdiag)
      mat=as.data.table(out$mat)
  }
  if (calc.differences==T) {
      cat("fast binless difference\n")
      ret=binless:::fast_binless_difference(out,1,alpha,lam2,lam1,tol_val,
                                            compute_patchnos=F,maxdiag=maxdiag)
      ret=as.data.table(ret)
  }
  mat[,patchno:=NULL]
  end.time = Sys.time()
  time.norm = diff(c(start.time, end.time))
  units(time.norm) = "mins"
  cat("*** Fast binless time=",time.norm,"mins\n")
  
  
  cat("postprocessing\n")
  start.time = Sys.time()
  #save(out, file=paste0(out_prefix,"_fast_out.RData"))
  #rm(out)
  #if (calc.differences==T) {
  #    save(ret, file=paste0(out_prefix,"_fast_ret.RData"))
  #}
  
  if (calc.differences==F) {
      lam1.sig=mat[distance>far.cutoff,median(log(signal))]
      if (is.na(lam1.sig)) {
        cat("Warning: using lam1.final as lambda1 estimate\n")
        lam1.sig=lam1.final
      }
      cat("*** Using signal lambda1=",lam1.sig,"\n")
      mat[,signal:=pmax(signal/exp(lam1.sig),1)]
  } else {
      lam1.sig=mat[distance>far.cutoff,median(log(signal)),by=name]$V1
      cat("*** Using signal lambda1=",lam1.sig,"\n")
      mat[name==name[1],signal:=pmax(signal/exp(lam1.sig[1]),1)]
      mat[name==name[.N],signal:=pmax(signal/exp(lam1.sig[2]),1)]
  
      # Thresholding is not obvious to compute, output both
      cat("*** Using difference lambda1=",lam1.final,"\n")
      ret[pos2-pos1>=maxdiag*base.res,difference:=1]
      ret[,signif:=log(difference)]
      ret[,signif:=ifelse(signif>0,pmax(signif-lam1.final, 0),
                                   pmin(signif+lam1.final, 0))]
      ret[,difference:=ifelse(difference>1,pmax(difference/exp(0.01), 1),
                          pmin(difference*exp(0.01), 1))]
  
  }
  
  end.time = Sys.time()
  time.save = diff(c(start.time, end.time))
  units(time.save) = "mins"
  cat("*** Postprocessing time=",time.save,"mins\n")
  
  
  cat("writing matrices\n")
  start.time = Sys.time()
  nbins = mat[,nlevels(bin1)] #number of bins
  rn = mat[,sort(as.numeric(tstrsplit(levels(bin1),"[[,]")[[2]]))]
  rn = paste0(rn,"-",rn+base.res)
  if (calc.differences==F) {
      # Save raw matrix
      newmat = mat[unclass(bin2)-unclass(bin1) <= maxdiag] #select diagonal
      newmat1 = sparseMatrix(i=newmat[,unclass(bin1)],j=newmat[,unclass(bin2)],
                             x=newmat[,observed], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_raw_band.RData")
      save(newmat1, file=outname)
  
      # Save binless matrix
      newmat1 = sparseMatrix(i=newmat[,unclass(bin1)],j=newmat[,unclass(bin2)],x=newmat[,binless],
                             dims=rep(nbins,2), dimnames=list(rn,rn), symmetric=T, index1=T)
  
      outname = paste0(out_prefix,"_binless_band.RData")
      save(newmat, file=outname)

      #Save signal matrix
      newmat = mat[(signal!=1)]
      newmat = sparseMatrix(i=newmat[,unclass(bin1)],j=newmat[,unclass(bin2)],x=newmat[,log10(signal)],
                         dims=rep(nbins,2), dimnames=list(rn,rn), symmetric=T, index1=T)
      outname = paste0(out_prefix,"_log10signal.RData")
      save(newmat, file=outname)
  } else {
      # Save raw matrix
      newmat = mat[unclass(bin2)-unclass(bin1) <= maxdiag] #select diagonal
      newmat1 = sparseMatrix(i=newmat[name==name[1],unclass(bin1)],j=newmat[name==name[1],unclass(bin2)],
                             x=newmat[name==name[1],observed], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_ref_raw_band.RData")
      save(newmat1, file=outname)
      #
      newmat2 = sparseMatrix(i=newmat[name==name[.N],unclass(bin1)],j=newmat[name==name[.N],unclass(bin2)],
                             x=newmat[name==name[.N],observed], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_target_raw_band.RData")
      save(newmat2, file=outname)
  
      # Save binless matrix
      newmat1 = sparseMatrix(i=newmat[name==name[1],unclass(bin1)],j=newmat[name==name[1],unclass(bin2)],
                             x=newmat[name==name[1],binless], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_ref_binless_band.RData")
      save(newmat1, file=outname)
      #
      newmat2 = sparseMatrix(i=newmat[name==name[.N],unclass(bin1)],j=newmat[name==name[.N],unclass(bin2)],
                             x=newmat[name==name[.N],binless], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_target_binless_band.RData")
      save(newmat2, file=outname)
  
      #save signal matrices
      newmat = mat[(signal!=1)]
      newmat1 = sparseMatrix(i=newmat[name==name[1],unclass(bin1)],j=newmat[name==name[1],unclass(bin2)],
                             x=newmat[name==name[1],log10(signal)], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_ref_log10signal.RData")
      save(newmat1, file=outname)
      #
      newmat2 = sparseMatrix(i=newmat[name==name[.N],unclass(bin1)],j=newmat[name==name[.N],unclass(bin2)],
                             x=newmat[name==name[.N],log10(signal)], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
  
      outname = paste0(out_prefix,"_target_log10signal.RData")
      save(newmat2, file=outname)
  
      #Save unthresholded difference matrix
      nbins = ret[,nlevels(bin1)] #number of bins
      rn = mat[,sort(as.numeric(tstrsplit(levels(bin1),"[[,]")[[2]]))]
      rn = paste0(rn,"-",rn+base.res)
      newmat = ret[difference!=1]
      newmat1 = sparseMatrix(i=newmat[name==name[.N],unclass(bin1)],j=newmat[name==name[.N],unclass(bin2)],
                             x=newmat[name==name[.N],log10(difference)], dims=rep(nbins,2), dimnames=list(rn,rn),
                             symmetric=T, index1=T)
      outname = paste0(out_prefix,"_log10diff.RData")
      save(newmat1, file=outname)
  
      #Save thresholded difference matrix
      newmat = ret[signif!=0]
      newmat2 = sparseMatrix(i=newmat[name==name[.N],unclass(bin1)],j=newmat[name==name[.N],unclass(bin2)],
                            x=newmat[name==name[.N],signif/log(10)], dims=rep(nbins,2), dimnames=list(rn,rn),
                            symmetric=T, index1=T)
      outname = paste0(out_prefix,"_log10diff_signif.RData")
      save(newmat2, file=outname)
  
  }
  
}

#' Read matrices written by \code{\link{chromosome_binless}}
#'
#' This procedure produces a mat object amenable to \code{\link{plot_binless_matrix}}.
#'
#' @param sample The name to use for the sample
#' @param sample2 If difference matrices are to be loaded, the name to use for
#'   the second sample
#' @param zoom_region If "full" (default), will load complete dataset. If a
#'   vector c(start,end) is passed, will read only this portion.
#' @param out_prefix The output prefix used to store matrices.
#' @param read.raw Whether to load _raw_band.RData file(s). Default TRUE.
#' @param read.binless Whether to load _binless_band.RData file(s). Default
#'   TRUE.
#' @param read.signal Whether to load _log10signal.RData file(s). Default TRUE.
#' @param read.difference Whether to load _log10diff.RData file (if applicable).
#'   Default TRUE.
#' @param read.signif Whether to load _log10diff_signif.RData file (if
#'   applicable). Default TRUE.
#'
#' @return mat data.table
#' @export
#'
#' @examples " "
read_chromosome_binless = function(sample, sample2=NA, zoom_region="full",
                              out_prefix="chromosome_binless", read.raw=T, read.binless=T,
                              read.signal=T, read.difference=T, read.signif=T) {
  calc.differences = !is.na(sample2)
  
  #prepare file loading
  sig.mask=c(read.raw,read.binless,read.signal)
  fnames=paste0(c("_raw_band","_binless_band","_log10signal"),".RData")[sig.mask]
  varnames=c("observed","binless","signal")[sig.mask]
  transnames=c("id","id","e10")[sig.mask]
  if (calc.differences==F){
    fnames=data.table(sample=sample,fname=paste0(out_prefix,fnames),variable=varnames,trans=transnames)
  } else {
    a=data.table(sample=sample,fname=paste0(out_prefix,"_ref",fnames),variable=varnames,trans=transnames)
    a=rbind(a,data.table(sample=sample2,fname=paste0(out_prefix,"_target",fnames),variable=varnames,trans=transnames))
    diff.mask=c(read.difference,read.signif)
    fnames=paste0(c("_log10diff","_log10diff_signif"),".RData")[diff.mask]
    varnames=c("difference","signif")[diff.mask]
    transnames=c("e10","e10")[diff.mask]
    fnames=rbind(a,data.table(sample=sample2,fname=paste0(out_prefix,fnames),variable=varnames,trans=transnames))
  }
  
  #read base resolution
  first_bin=rownames(get(load(fnames[1,fname])))[1]
  base.res=as.integer(tstrsplit(first_bin,"-"))
  base.res=base.res[2]-base.res[1]
  
  #load files and dcast
  mat = foreach(idx=1:(fnames[,.N]),.combine=rbind) %do% {
    mat=get(load(fnames[idx,fname]))
    if (nnzero(mat)>0) {
      rn=as.integer(tstrsplit(rownames(mat),"-")[[1]])
      mat=as.data.table(summary(mat))[,.(name=fnames[idx,sample],begin1=rn[as.integer(i)],begin2=rn[as.integer(j)],
                                         variable=fnames[idx,variable],value=x)]
      #
      if (!(length(zoom_region)==1 && zoom_region == "full")) {
        mat = mat[begin1>=zoom_region[1]&begin1<=zoom_region[2]&begin2>=zoom_region[1]&begin2<=zoom_region[2]]
      }
      #
      if (fnames[idx,trans] == "e10") {
        mat[,value:=10**value]
      } else {
        stopifnot(fnames[idx,trans]=="id")
      }
    } else {
      mat=data.table()
    }
    mat
  }
  mat = dcast(mat,name+begin1+begin2~variable,value.var="value")
  
  #polish columns and return
  mat[,name:=ordered(name,ifelse(calc.differences==T,c(sample,sample2),sample))]
  bin_labels=mat[,sort(unique(begin1))]
  bin_labels=paste0("[",bin_labels,",",bin_labels+base.res,")")
  mat[,bin1:=ordered(paste0("[",begin1,",",begin1+base.res,")"),bin_labels)]
  mat[,bin2:=ordered(paste0("[",begin2,",",begin2+base.res,")"),bin_labels)]
  setcolorder(mat,c("name","bin1","begin1","bin2","begin2"))
  mat
}
  



