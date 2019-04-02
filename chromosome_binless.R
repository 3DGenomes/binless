#!/usr/bin/env Rscript

library(binless)

args=commandArgs(trailingOnly=TRUE)
if (length(args)<= 1 || length(args) > 3) stop(paste0("\nargs: resolution sample1 [sample2]\n"))
resolution=as.integer(args[1])
sample=args[2]   #provide path to csdata object in RData format
if (length(args)==3) {
    sample2 = args[3] 
} else {
    sample2 = NA
}



### fixed parameters

ncores = 8
far.cutoff=10e6 #distance at which we should stop computing signal 
                #(used for final lambda1 estimation, to lower memory footprint,
                #and to write matrices)
out_prefix="chromosome_binless"
#zoom_region="full"
zoom_region=c(20000000,30000000)

chromosome_binless(sample, sample2=sample2, zoom_region=zoom_region, ncores=ncores, base.res=resolution,
                   far.cutoff=far.cutoff, out_prefix="chromosome_binless")

