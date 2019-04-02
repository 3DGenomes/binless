#!/usr/bin/env Rscript

library(doParallel)
library(data.table)
library(foreach)
library(binless)
library(Matrix)

makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

di_index = function(hicmat,nw,n,triag=T) {
  if(triag) {
    A = matrix(NA,n,n)
    A[lower.tri(A,diag=TRUE)] = hicmat
    A = t(A)
    A = makeSymm(A)
  } else {
    A = hicmat
  }
  #signal1 = rep(0,n)
  signal1 = foreach (i=(1+nw):(n-nw),.combine=c) %do% {
    vect_left = foreach (k=(i-nw):(i-1),.combine=c) %do% {
      if(k < 1) return(NA)
      if(A[i,k] > 0) return(log(A[i,k]))
      else return(0);
    }
    vect_left = vect_left[!is.na(vect_left)]
    vect_right = foreach (k=(i+1):(i+nw),.combine=c) %do% {
      if(k > n) return(NA)
      if(A[i,k] > 0) return(log(A[i,k]))
      else return(0);
    }
    vect_right = vect_right[!is.na(vect_right)]
    
    if(sum(vect_left) != 0 && sum(vect_right) != 0) return(t.test(vect_right,vect_left, paired = TRUE)$statistic)
    else return(0)
  }
  return(sd(signal1))
}


###variable arguments

args=commandArgs(trailingOnly=TRUE)
if (length(args)<= 1 || length(args) > 3) stop(paste0("\nargs: resolution sample1 [sample2]\n"))
resolution=as.integer(args[1])
sample=args[2]   #provide path to csdata object in RData format
if (length(args)==3) {
    calc.differences=T
    sample2 = args[3] 
} else {
    calc.differences=F
}



### fixed parameters

ncores = 16
nperf=100 #optimized binless perf iteration count
far.cutoff=10e6 #distance at which we should stop computing signal 
                #(used for final lambda1 estimation, to lower memory footprint,
                #and to write matrices)
out_prefix="chromosome_binless"
destdir="."
write.intermediate=F
write.logfile=F
zoom_region="full"
#zoom_region=c(20000000,40000000)
### deduced parameters
maxdiag = far.cutoff/resolution
submat_step = 20*resolution
submat_step_bin=floor(submat_step/resolution)
submat_length = 120*resolution
base.res=resolution #base resolution for the fused lasso signal detection


###remaining untouched parameters
if (exists("bf_per_kb")) bf_per_kb=as.integer(bf_per_kb) else bf_per_kb=50
if (exists("bf_per_decade")) bf_per_decade=as.integer(bf_per_decade) else bf_per_decade=10
if (exists("bins_per_bf")) bins_per_bf=as.integer(bins_per_bf) else bins_per_bf=10
if (exists("ngibbs")) ngibbs=as.integer(ngibbs) else ngibbs=35 #maximum number of iterations
if (exists("bg.steps")) bg.steps=as.integer(bg.steps) else bg.steps=5 #maximum number of steps where only the background model is fitted
if (exists("n_iter")) n_iter=as.integer(n_iter) else n_iter=100
if (exists("tol")) tol=as.numeric(tol) else tol=0.1 #relative tolerance on computed quantities upon convergence
if (exists("nperf")) nperf=as.integer(nperf) else nperf=75
if (exists("min.lambda2")) min.lambda2 = as.numeric(min.lambda2) else min.lambda2=0.1


if (exists("nouter")) nouter=as.integer(nouter) else nouter=25
if (exists("tol_val")) tol_val=as.integer(tol_val) else tol_val=0.2
if (exists("bg_steps")) bg_steps=as.integer(bg_steps) else bg_steps=5
if (exists("free_decay")) free_decay=as.integer(free_decay) else free_decay=10000

if (exists("ncores")) ncores=as.integer(ncores) else ncores=4 #parallelize on so many processors

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
    matsub = binless:::bin_data(csa,resolution=resolution)
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
    matsub = binless:::bin_data(csa,resolution=resolution)
    sddi=mean(matsub[,di_index(observed,floor(submat_length/resolution/10),floor((.N*2)**0.5)),by=name]$V1)
    st=cbind(counts_stats[i],data.table(std.directionality=sddi))
    st
  }
  di=di[order(-std.directionality)]
} else {
    #differences: compute sum of absolute differences larger than 0.1 on all submatrices
  di=foreach (i=seq_len(counts_stats[,.N]), csa=cs_arr, .errorhandling = 'remove',
              .packages=c("foreach","data.table","binless"), .combine=rbind) %dopar% {
     matsub = binless:::bin_data(csa,resolution=resolution)
     alpha=0.1
     lam2=2.5
     lam1=0
     out=binless:::fast_binless(matsub, matsub[,nlevels(bin1)],
                                alpha, lam2, lam1, nouter, tol_val,
                                bg_steps, free_decay, compute_patchnos = F,
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
di[,is.selected:=is.selected&cumsum(is.selected)<=10]
cs_arr = cs_arr[di[,is.selected==T]]
start_positions = start_positions[di[,is.selected==T]]


### run parallel normalizations and signal (difference) detections
### and extract sensible values for lambda1, lambda2 and alpha
all_chunks_norm = foreach (start_pos=start_positions, csb=cs_arr, .errorhandling = 'remove',
            .packages=c("foreach","data.table","binless"), .combine=rbind) %dopar% {
      matsub = binless:::bin_data(csb,resolution=resolution)
      #csb = cs_arr[[c]]
      csb <- normalize_binless(csb, ngibbs = ngibbs, ncores = 2, base.res = base.res, bg.steps = bg.steps, tol = tol,
                               bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=n_iter,
                               min.lambda2 = min.lambda2)
      save(csb,file=paste0(format(as.integer(start_pos),scientific = FALSE),"-",format(as.integer(start_pos+submat_length),scientific = FALSE),'.RData'))
      if (has_converged(csb, "signal")) {
        csb=bin_all_datasets(csb, ncores = 2)
        if (calc.differences==F) {
            csb=detect_binless_interactions(csb, ncores = 2, nperf = nperf)
        } else {
            csb=detect_binless_differences(csb, cs@experiments[1,name],
                                           ncores = 2, nperf = nperf)
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
saveRDS(all_chunks_norm,file=paste0(out_prefix,"_submatrix_parameters.rda"))
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
mat=binless:::bin_data(cs,resolution=resolution)


cat("fast binless\n")

if (write.intermediate==T) {
    csvdata=paste0(out_prefix,"_tmp.csv.gz")
    out=binless:::fast_binless(mat, mat[,nlevels(bin1)], alpha, lam2, lam1,
                               nouter, tol_val, bg_steps, free_decay, compute_patchnos = F,
                               csv_out=csvdata, maxdiag=maxdiag)
    df=read.csv(csvdata)
    mat = mat[,c('name','bin1','pos1','bin2','pos2','distance','observed','nobs')]
    mat = cbind(mat[,c('name','bin1','pos1','bin2','pos2','distance','observed','nobs')],df)
    out$mat=mat
    rm(df)
    unlink(csvdata)
} else {
    out=binless:::fast_binless(mat, mat[,nlevels(bin1)], alpha, lam2, lam1,
                               nouter, tol_val, bg_steps, free_decay, compute_patchnos = F,
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
save(out, file=paste0(out_prefix,"_fast_out.RData"))
rm(out)
if (calc.differences==T) {
    save(ret, file=paste0(out_prefix,"_fast_ret.RData"))
}

if (calc.differences==F) {
    lam1.sig=mat[distance>far.cutoff,median(log(signal))]
    cat("*** Using signal lambda1=",lam1.sig,"\n")
    mat[,signal:=pmax(signal/exp(lam1.sig),1)]
} else {
    lam1.sig=mat[distance>far.cutoff,median(log(signal)),by=name]$V1
    cat("*** Using signal lambda1=",lam1.sig,"\n")
    mat[name==name[1],signal:=pmax(signal/exp(lam1.sig[1]),1)]
    mat[name==name[.N],signal:=pmax(signal/exp(lam1.sig[2]),1)]

    # Thresholding is not obvious to compute, output both
    cat("*** Using difference lambda1=",lam1.final,"\n")
    ret[pos2-pos1>=maxdiag*resolution,difference:=1]
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
rn = mat[,paste0(sort(as.numeric(tstrsplit(levels(bin1),"[[,]")[[2]])))] #row names in quique format
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
    rn = ret[,paste0(sort(as.numeric(tstrsplit(levels(bin1),"[[,]")[[2]])))] #row names in quique format
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



