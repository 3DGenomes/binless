#' @include binless.R
NULL

#' Apply the ICE algorithm to a binned matrix
#' 
#' @param mat a matrix obtained by grouping
#' @param niterations positive integer. Number of iterations to perform
#'
#' @return a CSbinned object containing the ICEd matrix
#' @export
#'
#' @examples
iterative_normalization = function(mat, niterations=100, namecol="name", verbose=T) {
  if (verbose==T) cat("*** iterative normalization with ",niterations," iterations\n")
  raw=mat[,.(name,bin1,bin2,observed)]
  setkey(raw,name,bin1,bin2)
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
  if ("begin1" %in% names(mat)) 
    binned=merge(binned,mat[,.(name,bin1,begin1,end1,bin2,begin2,end2)],
                 by=c("name","bin1","bin2"),all.x=T)
  setkey(binned,name,bin1,bin2)
  binned
}

