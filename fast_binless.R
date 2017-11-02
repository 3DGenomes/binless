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


#In fast binless, you must ensure no counter diagonal nor row/column is completely zero
#if there are not too many, you can add 1 to the observed counts
#also, because of the simplified decay, the very last bin cannot be zero in all datasets
#empty rows:
rbind(mat[,.(bin=bin1,observed)],mat[,.(bin=bin2,observed)])[,.(sum(observed)),by=bin][V1==0]
#empty counter diagonals:
mat[,.(d=unclass(bin2)-unclass(bin1),observed)][,sum(observed),by=d][V1==0]
#last bin:
mat[unclass(bin1)==1&unclass(bin2)==max(unclass(bin2))]

#The fast binless algorithm computes a binless normalization without estimating
#the fusion penalty and the significance threshold.
#We do a maximum of nouter steps (or less if the relative precision is lower than tol_val)
#and hold the fusion penalty fixed at lam2.
#Play with lam2 to see its effect. This is the parameter that is optimized in the
#full-blown binless, along with another threshold (lambda1, set to zero here),
#which determines the significance of a given signal/difference contribution
nouter=25
lam2=5
tol_val=2e-1
bg_steps=5
out=binless:::fast_binless(mat, mat[,nlevels(bin1)], lam2, nouter, tol_val, bg_steps)


#Here follow the plots of the observed and fitted quantities (be sure to check out signal and binless plots)
#all data
a=as.data.table(out$mat)
#observed matrix (input data)
plot_binless_matrix(a, upper="observed", lower="observed")
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_fast_binless_observed.pdf", width=18,height=8)
#fitted background
plot_binless_matrix(a, upper="log_background", lower="log_background", trans="identity")
#fitted biases
ggplot(data.table(bin=1:nlevels(mat[,bin1]),log_biases=out$log_biases))+geom_point(aes(bin,log_biases,colour="cpp"))
#biases matrix
plot_binless_matrix(a, upper="log_biases", lower="log_biases", trans="identity")
#fitted decay
ggplot(data.table(distance=1:nlevels(mat[,bin1]),log_decay=out$log_decay))+geom_point(aes(distance,log_decay,colour="cpp"))
#decay matrix
plot_binless_matrix(a, upper="log_decay", lower="log_decay", trans="identity")
#signal matrix ( = what is different from the background)
plot_binless_matrix(a, upper="log_signal", lower="log_signal", trans="identity")
#weights ( = 1/variance )
plot_binless_matrix(a, upper="weights", lower="weights")
#binless matrix ( = signal + decay)
plot_binless_matrix(a, upper="log_binless", lower="log_binless", trans="identity")
#binless and observed
plot_binless_matrix(a, upper="log_binless", lower="log(observed)", trans="identity")
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_fast_binless.pdf", width=18,height=8)


### Difference detection

#now we compute differences between the two datasets
ref=1
lam2=5
tol_val=2e-1
diff=as.data.table(binless:::fast_binless_difference(out, lam2, ref, tol_val))

#log(observed)
plot_binless_matrix(diff, upper="observed", lower="observed")
#log difference of all datasets wrt ref
plot_binless_matrix(diff[name!=ref], upper="log_difference/log(2)", lower="log_difference/log(2)", trans="identity", label="log2 FC", limits = c(-3,3))
#ggsave(filename="example/rao_HiCall_FOXP1ext_2.3M_fast_binless_difference.pdf", width=10,height=8)

