#' binless: A package for resolution-independent Hi-C normalization
#'
#' This package performs binless normalization, a negative binomial generalized
#' additive model regression that uses dangling and rejoined ends to normalize
#' Hi-C data in a resolution-independent way. It also allows for interaction
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
#' \code{\link{normalize_binless}}
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
#' @useDynLib binless, .registration = TRUE
#' @import methods
#' @import data.table
#' @import doParallel
#' @import foreach
#' @import matrixStats
#' @import ggplot2
#' @import MASS
#' @import Matrix
#' @import quadprog
#' @importFrom scales muted
#' @importFrom Hmisc cut2
#' @importFrom dplyr ntile
#' @importFrom utils packageVersion
#'
#' @docType package
#' @name binless
NULL

