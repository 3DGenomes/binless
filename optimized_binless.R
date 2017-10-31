library(ggplot2)
library(scales)
library(data.table)
library(binless)

### Normalization using binless, optimizing the shrinkage and thresholding penalties

#Data must have been preprocessed, and a CSnorm object must be created and stored in a variable called cs.
#load("example/rao_HiCall_FOXP1ext_csnorm.RData")


### Run the normalization

#default parameters should be correct for 6-cutter datasets
ncores=4 #parallelize on so many processors
ngibbs=15 #maximum number of iterations
base.res=5000 #base resolution for the fused lasso signal detection
bg.steps=5 #maximum number of steps where only the background model is fitted
tol=1e-1 #relative tolerance on computed quantities upon convergence
cs <- normalize_binless(cs, ngibbs = ngibbs, ncores = ncores, base.res = base.res, bg.steps = bg.steps, tol = tol)

# If you get the message "Normalization has converged", congratulations! You can save the normalization like so
#save(cs,file="example/rao_HiCall_FOXP1ext_csnorm_optimized.RData")
# Alternately, if you want to save space and drop all intermediates that can be recomputed,
# and all diagnostics information (keeping just the normalization parameters) write
#save_stripped(cs, file="example/rao_HiCall_FOXP1ext_csnorm_optimized_stripped.RData")
# You can jump to the next section.

#Otherwise, we will first see how the normalization went so far
#We first plot the four log-likelihoods / cross-validation scores, which should reach a plateau
plot_diagnostics(cs)$plot 
#We can also see how all parameters behave
plot_diagnostics(cs)$plot2

#If you see damped oscillations, most likely you only need to extend the normalization by a couple of steps, say 10
cs <- normalize_binless(cs, restart=T, ngibbs = 10, ncores = ncores)

#Otherwise, try the following
# - change the base resolution. Watch out: the algorithm CPU and memory usage scale quadratically with the base resolution
# - normalize a larger portion of the genome. Go through all preprocessing but extend your region of interest by, say, 1Mb on each side.


### Obtaining binned matrices

#If you do not have a normalized dataset already, load it from the disk
#load("example/rao_HiCall_FOXP1ext_csnorm_optimized.RData")
#load_stripped("example/rao_HiCall_FOXP1ext_csnorm_optimized_stripped.RData", ncores=ncores)

#first, we create all binned matrices, and prepare all intermediates for binless calculation
#in each of the following calls, there is a resolution parameter, which you can set to
#something different than base.res. However, computing binless matrices at a different
#resolution than the one used during normalization is not advised.
cs=bin_all_datasets(cs, ncores=ncores)

#You can get the generated matrices (as data.table objects) with the call
mat=get_binned_matrices(cs)
#Here we plot the observed counts
plot_binless_matrix(mat, upper="log(observed)", lower="log(observed)")
#The optimized background (upper panel), and its estimated standard deviation 
plot_binless_matrix(mat, upper="log(background)", lower="log(background.sd)")
#The "normalized" (ICE-like) matrix, where we simply subtract the estimated biases from the raw data 
plot_binless_matrix(mat, upper="log(normalized)", lower="log(normalized)")

#If you have more than one dataset, you can group data, by providing the appropriate grouping
#Imagine you have two replicates and two conditions, you could group the data using group="condition"
#to generate one combined matrix for each condition. Note that the previous call to bin_all_datasets
#is equivalent to calling group_datasets with group="all"
cs=group_datasets(cs, group="condition", ncores=ncores)


### Obtaining binless signal matrices

#We perform the final computation of the binless signal matrices for each dataset (or group of datasets)
cs=detect_binless_interactions(cs, ncores=ncores)

#You can get the generated matrices (as data.table objects)
#For commodity, binned matrices are reported as well
mat=get_binless_interactions(cs)
#We plot the fused lasso signal that is estimated to be significantly different from the background model
#"signal" is a fold change
plot_binless_matrix(mat, upper="log2(signal)", lower="log2(signal)")
#The binless matrix corresponds to the signal matrix plus the estimated decay
plot_binless_matrix(mat, upper="log(binless)", lower="log(binless)")
#Binless matrices can be visually compared to the raw data to validate the normalization
plot_binless_matrix(mat, upper="log(binless)", lower="log(observed)")


### Obtaining binless difference matrices

#To compute differences between datasets (or groups of datasets) you must choose a reference
#Have a look at the cs@experiments data.table. If you did not perform any grouping take one of the
#entries in the "name" column. Otherwise, refer to that group's name
ref = cs@experiments[1,name]

#We then perform the difference computation for each remaining dataset (or group of datasets)
cs=detect_binless_differences(cs, ref = ref, ncores=ncores)

#You can get the generated matrices (as data.table objects)
#For commodity, binned matrices are reported as well
mat=get_binless_differences(cs, ref = ref)
#We plot the difference that is estimated to be significant between each dataset and the reference
#log difference of all datasets wrt ref
#"difference" is a fold change
plot_binless_matrix(mat[name!=ref], upper="log2(difference)", lower="log2(difference)")

### This concludes the quick tour on binless normalization!
### We hope you will enjoy the software and please do not hesitate to report
### any feature, bug or question to our GitHub issue tracker

