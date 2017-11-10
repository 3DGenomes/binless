% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_binless.R
\name{detect_binless_interactions}
\alias{detect_binless_interactions}
\title{Perform binless interaction detection using fused lasso}
\usage{
detect_binless_interactions(cs, resolution = cs@settings$base.res,
  group = "all", ncores = 1, tol.val = cs@settings$tol, verbose = T,
  fix.lambda1 = F, fix.lambda1.at = NA, fix.lambda2 = F,
  fix.lambda2.at = NA)
}
\arguments{
\item{cs}{}

\item{resolution}{}

\item{group}{}

\item{ncores}{}

\item{ref}{}

\item{niter}{number of IRLS iterations, and BIC iterations within}
}
\description{
Perform binless interaction detection using fused lasso
}