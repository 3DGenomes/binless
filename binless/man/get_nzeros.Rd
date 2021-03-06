% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{get_nzeros}
\alias{get_nzeros}
\title{count number of zeros in a given cut site, distance bin and signal bin}
\usage{
get_nzeros(cs, sbins, ncores = 1)
}
\value{
a data table, keyed on name,dbin,id1,pos1,bin1,bin2,dir,cat with the following columns
- name: name of the dataset
- id1, pos1, bin1: coordinates of the cut site
- bin2, dbin: signal/distance bin in which we looked at the intersections
- dir: fwd (rev) for contacts with downstream (resp upstream) cut-sites
- cat: whether we consider contacts on the left (contact L) or on the right (contact R) of this cut site
- ncross: number of cut site intersections (crossings) in this signal/distance/direction/category bin.
 We discard anything below cs@settings$dmin. Note that sum(ncross) is four times the total number of crossings (per dataset).
- nnz: number of non-zero contacts in this signal/distance/direction bin (max 2 per crossing).
 Note that sum(nnz) is twice the number of nonzeros (per dataset)
- nzero: number of zeros in this signal/distance/direction bin. We have nzero = 2*ncross - nnz.
 Note that sum(nzero) is twice the number of zeros (per dataset)
 For speed purposes downstream, we only return the entries where nzero>0, which corresponds to most of the entries anyway.
 We therefore have sum(nnz+nzero) equal to approximately twice the number of detectable counts, e.g. 8x the number of crossings
}
\description{
count number of zeros in a given cut site, distance bin and signal bin
}
\keyword{internal}
