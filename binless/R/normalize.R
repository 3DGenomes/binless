#' @include binless.R
NULL

#' Convert a sparse counts data.table to a dense one by adding rows with zero counts
#'
#' @param counts,biases data.tables as returned by \code{\link{prepare_for_sparse_cs_norm}}
#' @param biases2 data.table of biases for id2 column of counts. If NULL (default), use that of biases
#' @param circularize integer. Length of the genome if circular
#' @param dmin numeric. Minimum distance to be considered a contact
#'
#' @return a counts data.table with zeros filled according to cut sites provided in biases (and biases2 if available)
#' @keywords internal
#' @export
#' @section Warning:
#' Memory-intensive
#' @examples
fill_zeros = function(counts,biases,biases2=NULL,circularize=-1L,dmin=0) {
  if (is.null(biases2)) biases2=biases
  runname=biases[,unique(name)]
  if (length(runname)>1) {
    foreach (i=runname, .combine=rbind) %do%
      fill_zeros(counts[name==i],biases[name==i],biases2[name==i],circularize=circularize,dmin=dmin)
  } else {
    if (biases[,.N]==0 | biases2[,.N]==0) return(data.table())
    newcounts=CJ(biases[,paste(id,pos)],biases2[,paste(id,pos)])
    newcounts[,c("id1","pos1"):=tstrsplit(V1, " ")]
    newcounts[,c("id2","pos2"):=tstrsplit(V2, " ")]
    newcounts[,c("id1","id2","pos1","pos2","V1","V2"):=
                list(as.integer(id1),as.integer(id2),as.integer(pos1),as.integer(pos2),NULL,NULL)]
    newcounts=newcounts[pos1<pos2]
    setkey(newcounts, id1, id2, pos1, pos2)
    setkey(counts, id1, id2, pos1, pos2)
    newcounts=counts[newcounts]
    newcounts[is.na(contact.close),contact.close:=0]
    newcounts[is.na(contact.far),contact.far:=0]
    newcounts[is.na(contact.up),contact.up:=0]
    newcounts[is.na(contact.down),contact.down:=0]
    newcounts[is.na(name),name:=runname]
    if (circularize>0) {
      newcounts[,distance:=pmin(abs(pos2-pos1), circularize+1-abs(pos2-pos1))]
    } else {
      newcounts[,distance:=abs(pos2-pos1)]
    }
    newcounts=newcounts[distance>=dmin]
    newcounts
  }
}

#' Initial guess for normalization
#' @keywords internal
#' @export
#' 
optimize_stan_model = function(model, data, iter, verbose, init, ...) {
  out=capture.output(op<-optimizing(model, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=init, ...))
  cat(out,sep="\n")
  if (length(grep("Line search failed",tail(out,1)))>0) {
    op=optimizing(model, data=data, as_vector=F, hessian=F, iter=2, verbose=verbose, init=init, algorithm="Newton", ...)
    cat("!!! Line search error occurred, performed 2 steps of Newton optimization\n")
  }
  return(op)
}
