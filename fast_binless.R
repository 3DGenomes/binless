library(ggplot2)
library(scales)
library(data.table)
library(binless)

### Normalization using the fast and approximate binless algorithm

#Fast binless requires a binned matrix of raw data at high resolution.
#You can either load it yourself using data.table::fread (file format description in README.md)
#don't forget to use stringsAsFactors=T
#On a mac, use gzcat instead of zcat (or feed it the uncompressed file)
#mat=fread("zcat example/rao_HiCall_FOXP1ext_2.3M_mat_5kb.dat.gz",stringsAsFactors = T)

#You can also extract it from the CSnorm object built during preprocessing
#load("example/rao_HiCall_FOXP1ext_csnorm.RData")
mat=binless:::bin_data(cs,resolution=5000)
#write.table(mat,file = "example/rao_HiCall_FOXP1ext_2.3M_mat_5kb.dat", quote = T, row.names = F)


#The fast binless algorithm computes a binless normalization without estimating
#the fusion penalty and the significance threshold.
#We do a maximum of nouter steps (or less if the relative precision is lower than tol_val)
#and hold the fusion penalty fixed at lam2.
#Play with lam2 to see its effect. This is the parameter that is optimized in the
#full-blown binless, along with another threshold (lambda1, set to zero here),
#which determines the significance of a given signal/difference contribution
nouter=25
lam2=5
alpha=10
tol_val=2e-1
bg_steps=5
out=binless:::fast_binless(mat, mat[,nlevels(bin1)], lam2, alpha, nouter, tol_val, bg_steps)


#Here follow the plots of the observed and fitted quantities (be sure to check out signal and binless plots)
#all data
a=as.data.table(out$mat)
#observed matrix (input data)
plot_binless_matrix(a, upper="observed", lower="observed")
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_fast_binless_observed.pdf", width=18,height=8)
#number of observables (input data)
plot_binless_matrix(a, upper="nobs", lower="nobs")
#fitted background
plot_binless_matrix(a, upper="background", lower="background")
#fitted biases
ggplot(data.table(bin=1:nlevels(mat[,bin1]),log_biases=out$log_biases))+geom_point(aes(bin,log_biases))
#biases matrix
plot_binless_matrix(a, upper="biasmat", lower="biasmat")
#fitted decay
ggplot(unique(a[,.(distance,decaymat)]))+geom_line(aes(distance,decaymat))+scale_x_log10()+scale_y_log10()
#decay matrix
plot_binless_matrix(a, upper="decaymat", lower="decaymat")
#signal matrix ( = what is different from the background)
plot_binless_signal_matrix(a)
#weights ( = 1/variance )
plot_binless_matrix(a, upper="weight", lower="weight")
#binless matrix ( = signal + decay)
plot_binless_matrix(a, upper="binless", lower="binless")
#binless and observed
plot_binless_matrix(a, upper="binless", lower="observed")
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_fast_binless.pdf", width=18,height=8)


### Difference detection

#now we compute differences between the two datasets
ref=1
lam2=5
alpha=1
tol_val=2e-1
diff=as.data.table(binless:::fast_binless_difference(out, lam2, ref, alpha, tol_val))

#log(observed)
plot_binless_matrix(diff, upper="observed", lower="observed")
#log difference of all datasets wrt ref
plot_binless_difference_matrix(diff[name!=ref])
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_fast_binless_difference.pdf", width=10,height=8)

