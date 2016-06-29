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
                       data="data.table",
                       biases="data.table",
                       counts="data.table"))

#' Class to hold cut-site normalization data
#'
#' @slot experiments data.table. 
#' @slot design data.table. 
#' @slot biases data.table. 
#' @slot counts data.table. 
#' @slot par list. 
#' @slot diagnostics list. 
#' @slot pred list. 
#' @slot binned list. 
#'
#' @return
#' @export
#'
#' @examples
setClass("CSnorm",
         slots = list(experiments="data.table",
                      design="data.table",
                      biases="data.table",
                      counts="data.table",
                      par="list",
                      diagnostics="list",
                      pred="list",
                      binned="list"))
