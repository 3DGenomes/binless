library(csnorm)

setwd("/home/yannick/simulations/cs_norm")

args=commandArgs(trailingOnly = T)
prefix=args[1]
coverage=as.integer(args[2])
square.size=as.integer(args[3])
bf_per_kb=as.numeric(args[4])

message("normalization on ",prefix, "with coverage ",coverage, " square size ",square.size," and bf_per_kb ",bf_per_kb)

biases=fread(paste0("data/",prefix,"_biases.dat"))
setkey(biases,id)
counts=fread(paste0("data/",prefix,"_counts.dat"))
load(paste0("data/",prefix,"_meanfield_100.RData"))

oppar=run_split_parallel(counts, biases, square.size=square.size, coverage=coverage, bf_per_kb=bf_per_kb,
                         bf_per_decade=5, distance_bins_per_decade=100, verbose = T, iter=100000, ncpus=30, homogenize=F,
                         outprefix=paste0("tmp/",prefix,"_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb))
                         #circularize=4042929)
save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))
oppar=postprocess(biases, counts, meanfield, oppar, resolution=10000, ncores=30, predict.all.means=T)
oppar$ice=iterative_normalization(oppar$mat, niterations=1, resolution=10000, return.binned=T)
save(oppar, file = paste0("data/",prefix,"_op_maxcount_-1_parallel_inhomogeneous_cov",coverage,"X_sq",round(square.size/1000),"k_bfpkb",bf_per_kb,".RData"))


