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
#' \code{\link{run_gauss}}
#' \code{\link{run_exact}}
#'
#' @section Postprocessing:
#' \code{\link{bin_all_datasets}}
#' \code{\link{group_datasets}}
#' \code{\link{detect_binned_interactions}}
#' \code{\link{detect_binned_differences}}
#' \code{\link{detect_binless_interactions}}
#' \code{\link{detect_binless_differences}}
#' \code{\link{generate_genomic_biases}}
#' \code{\link{plot_binless_matrix}}
#'
#' @useDynLib csnorm, .registration = TRUE
#' @import methods
#' @import splines
#' @import data.table
#' @import doParallel
#' @import foreach
#' @import matrixStats
#' @import ggplot2
#' @import MASS
#' @import Matrix
#' @import quadprog
#' @import igraph
#' @importFrom Hmisc cut2
#' @importFrom dplyr ntile
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
  cat(" from ", nreads, " count reads\n", sep="")
  cat(" Reads density excl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  nreads=nreads+object@biases[,sum(dangling.L+dangling.R+rejoined)]
  cat(" Reads density incl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  if (object@data[,.N]>0) {
    cat("Original data has ",object@data[,.N], " reads categorized as follows\n", sep="")
    show(object@data[,.N,keyby=category])
  } else {
    cat("Original data not stored in object")
  }
})

#' Virtual class for one interaction detection
#' @keywords internal
setClass("CSinter",
         slots = list(mat="data.table",
                      par="list",
                      settings="list"),
         contains = "VIRTUAL")

#' Binned interaction detection
#' @keywords internal
setClass("CSsig", slot = list(threshold="numeric"), contains = "CSinter")
setMethod("show", signature="CSsig", definition=function(object) {
  cat("        Significant interactions wrt expected") 
  cat(" (threshold=", object@threshold,")\n")
})

#' Binless difference detection
#' @keywords internal
setClass("CSdiff", slots = list(threshold="numeric", ref="character"), contains = "CSinter")
setMethod("show", signature="CSdiff", definition=function(object) {
  cat("        Significant differences wrt ", object@ref)
  cat(" (threshold=", object@threshold,")\n")
})

#' Binless interaction detection
#' @keywords internal
setClass("CSbsig", slots = list(cts="data.table",state="list"), contains = "CSinter")
setMethod("show", signature="CSbsig", definition=function(object) {
    cat("        Binless interactions wrt expected\n") 
})

#' Binless difference detection
#' @keywords internal
setClass("CSbdiff",
         slots = list(ref="character",cts="data.table",cts.ref="data.table",state="list"),
         contains = "CSinter")
setMethod("show", signature="CSbdiff", definition=function(object) {
  cat("        Binless differences wrt ", object@ref, "\n") 
})

#' Class for one dataset grouping at a given (base) resolution
#'
#' @slot mat data.table. 
#' @slot interactions list. 
#' @slot resolution in bases.
#' @slot group 
#' @slot cts data.table. The counts and predicted means used in all calculations.
#' @slot par
#' @slot names 
#'
#' @return
#' @keywords internal
#' @export
#'
#' @examples
setClass("CSgroup",
         slots = list(mat="data.table",
                      interactions="list",
                      resolution="numeric",
                      group="character",
                      cts="data.table",
                      par="list",
                      names="data.table"))

setMethod("show",signature="CSgroup",definition=function(object) {
  if (length(object@group)==1 && object@group=="all") {
    cat("   *** Individual at", object@resolution/1000, "kb resolution") 
  } else {
    cat("   *** Group [", object@group,"] at", object@resolution/1000, "kb resolution")
  }
  cat( " : ", object@names[,do.call(paste,c(as.list(groupname), sep=" / "))], "\n")
  if (length(object@interactions)==0) {
    cat("        No interactions computed\n")
  } else {
    lapply(object@interactions, show)
  }
  cat("\n")
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
                      zeros="data.table",
                      diagnostics="list",
                      pred="list",
                      groups="list"))

setMethod("show",signature="CSnorm",definition=function(object) {
  cat("Cut-site normalization object\n")
  cat(" Number of datasets: ", object@experiments[,.N], "\n", sep="")
  if(object@settings$circularize>0) {
    cat(" Genome is circular with size ", object@settings$circularize, "\n", sep="")
  } else {
    cat(" Genome is linear\n", sep="")
  }
  cat(" There are ", object@biases[,uniqueN(pos)], " unique cut sites\n", sep="")
  ncounts=object@counts[,sum((contact.close>0)+(contact.far>0)+(contact.up>0)+(contact.down>0))]
  cat(" and ", ncounts, " nonzero counts\n", sep="")
  nreads=object@counts[,sum(contact.close+contact.far+contact.up+contact.down)]
  cat(" from ", nreads, " reads\n", sep="")
  cat(" First cut site at ", object@biases[,min(pos)], "\n", sep="")
  cat(" Last cut site at ", object@biases[,max(pos)], "\n", sep="")
  cat(" Cut sites span ", object@biases[,(max(pos)-min(pos))/1000], " kb\n", sep="")
  cat(" Reads density excl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  nreads=nreads+object@biases[,sum(dangling.L+dangling.R+rejoined)]
  cat(" Reads density incl. biases: ", round(nreads/object@biases[,max(pos)-min(pos)]*1000), " reads per kilobase (rpkb)\n", sep="")
  cat(" Experimental design matrix:\n")
  show(merge(object@experiments[,.(name,condition,replicate,enzyme)],object@design,by="name"))
  if (length(object@par)==0) {
    cat(" Dataset not yet normalized\n")
  } else {
    cat(" Normalized dataset\n")
    cat("  lambda_iota: ",object@par$lambda_iota, "\n  lambda_rho: ",object@par$lambda_rho, "\n  lambda_diag: ",object@par$lambda_diag,"\n")
    cat("  dispersion: ",object@par$alpha,"\n  log likelihood: ", object@par$value, "\n")
    if (cs@diagnostics$params[,.N]>0) {
      if (has_converged(cs)==T) cat("  Normalization has converged\n") else  cat("  WARNING: Normalization did not converge!\n")
    }
    #
    ngroups=length(object@groups)
    if (ngroups==0) {
      cat(" No groupings available")
    } else {
      cat("\n ### ", ngroups, " groupings available\n\n", sep="")
      lapply(object@groups, show)
    }
  }
})
