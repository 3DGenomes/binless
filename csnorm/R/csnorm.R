#' csnorm: A package for cut-site normalization
#'
#' This package performs cut-site normalization, a negative binomial generalized
#' additive model regression that uses dangling and rejoined ends to normalize
#' 3C-like data in a resolution-independent way. It also allows for interaction
#' or interaction difference detection. See publication TODO.
#' 
#' @section Preprocessing:
#' \code{\link{read_and_prepare}}
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
#'
#' @docType package
#' @name csnorm
NULL
