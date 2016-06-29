#' csnorm: A package for cut-site normalization
#'
#' This package performs cut-site normalization, a negative binomial generalized
#' additive model regression that uses dangling and rejoined ends to normalize
#' 3C-like data in a resolution-independent way. It also allows for interaction
#' or interaction difference detection. See publication TODO.
#' 
#' @section Preprocessing:
#' \code{\link{read_and_prepare}} high-level wrapper for preprocessing steps
#' \code{\link{read_tsv}} read TadBit .tsv file
#' \code{\link{categorize_by_new_type}} add category label based on cut-site information
#' \code{\link{prepare_for_sparse_cs_norm}} summarize reads information per cut-site
#' \code{\link{dset_statistics}} print general statistics on this dataset
#' \code{\link{generate_fake_dataset}} follows model specification
#'
#' @section Normalization:
#' \code{\link{run_split_parallel}}
#'
#' @section Postprocessing:
#' \code{\link{postprocess}}
#' \code{\link{iterative_normalization}}
#' \code{\link{thresholds_estimator}}
#'
#' @useDynLib csnorm, .registration = TRUE
#' @import rstan
#' @import data.table
#' @import Hmisc
#' @import doParallel
#' @import foreach
#' @import MASS
#' @import matrixStats
#' @import ggplot2
#' @import methods
#'
#' @docType package
#' @name csnorm
NULL


#' Class to hold a single experiment
#'
#' @slot info list. Information on this experiment
#' @slot settings list. Settings for the normalization procedures
#' @slot data data.table. The reads data for that experiment
#' @slot biases data.table. The cut-site information for that experiment
#' @slot counts data.table. The reads data for that experiment
#'
#' @return
#' @export
#'
#' @examples
setClass("CSdata",
          slots = list(info="list",
                       settings="list",
                       data="data.table",
                       biases="data.table",
                       counts="data.table"))

setMethod("show",signature="CSdata",definition=function(object) {
  cat("A dataset for cut-site normalization\n")
  cat(" Input file: ", object@info$filename, "\n", sep="")
  cat(" Dataset name: ", object@info$name, "\n", sep="")
  cat(" Experiment: ", object@info$experiment, "\n", sep="")
  cat(" Condition: ", object@info$condition, "\n", sep="")
  cat(" Replicate: ", object@info$replicate, "\n", sep="")
  cat(" Enzyme: ", object@info$enzyme, "\n", sep="")
  cat(" Dangling end positions: L: ", object@info$dangling.L, " R: ", object@info$dangling.R, "\n", sep="")
  cat(" Maxlen: ", object@info$maxlen, "\n", sep="")
  if(object@settings$circularize>0) {
    cat(" Genome is circular with size ", object@settings$circularize, "\n", sep="")
  } else {
    cat(" Genome is linear\n", sep="")
  }
  cat(" There are ", object@biases[,.N], " cut sites\n", sep="")
  ncounts=object@counts[,sum((contact.close>0)+(contact.far>0)+(contact.up>0)+(contact.down>0))]
  cat(" and ", ncounts, " nonzero counts\n", sep="")
  nreads=object@counts[,sum(contact.close+contact.far+contact.up+contact.down)]
  cat(" from ", nreads, " reads\n", sep="")
  cat(" Reads density excl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  nreads=nreads+object@biases[,sum(dangling.L+dangling.R+rejoined)]
  cat(" Reads density incl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  if (object@data[,.N]>0) show(object@data[,.N,keyby=category])
})

#' Class to hold cut-site normalization data
#'
#' @slot experiments data.table. Essential information for each experiment
#' @slot design data.table. A design matrix for the normalization
#' @slot settings list. Additional parameters for the normalization
#' @slot biases data.table. Merged bias information for all experiments
#' @slot counts data.table. Merged count information for all experiments
#' @slot par list. Parameters of normalization
#' @slot diagnostics list. Diagnostics output of parallel runs
#' @slot pred list. Predicted quantities
#' @slot binned list. Binned matrices and differences
#'
#' @return
#' @export
#'
#' @examples
setClass("CSnorm",
         slots = list(experiments="data.table",
                      design="data.table",
                      settings="list",
                      biases="data.table",
                      counts="data.table",
                      par="list",
                      diagnostics="list",
                      pred="list",
                      binned="list"))

setMethod("show",signature="CSnorm",definition=function(object) {
  cat("Cut-site normalization object\n")
  cat(" Number of datasets: ", object@experiments[,.N], "\n", sep="")
  if(object@settings$circularize>0) {
    cat(" Genome is circular with size ", object@settings$circularize, "\n", sep="")
  } else {
    cat(" Genome is linear\n", sep="")
  }
  cat(" There are ", object@biases[,.N], " cut sites\n", sep="")
  ncounts=object@counts[,sum((contact.close>0)+(contact.far>0)+(contact.up>0)+(contact.down>0))]
  cat(" and ", ncounts, " nonzero counts\n", sep="")
  nreads=object@counts[,sum(contact.close+contact.far+contact.up+contact.down)]
  cat(" from ", nreads, " reads\n", sep="")
  cat(" Reads density excl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  nreads=nreads+object@biases[,sum(dangling.L+dangling.R+rejoined)]
  cat(" Reads density incl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
})
