library(csnorm)
library(data.table)
library(flsa)

setwd("/home/yannick/simulations/cs_norm")
load("flsa_failure_data.RData")

ret = flsa(submat[,value], connListObj=cmat, verbose=T, splitCheckSize = 1e5)
save(ret, file="flsa_failure_ret.RData")

