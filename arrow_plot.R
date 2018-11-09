library(ggplot2)
library(scales)
library(data.table)
library(foreach)
library(binless)

### Arrow plot

#Here, we will show how to quickly generate base-resolution (arrow) plots.
#If you intend to do binless normalization, you should rather refer to
#preprocessing.R, which covers these plots as well.
#All functions described here are documented. If you use rstudio, simply
#hit ?function_name in the console to get its documentation.

#we read the input file and store it in a data.table
#if you use a Mac, use gzcat instead of zcat, or provide the path to an uncompressed file
#refer to README.md for a description of the tsv file format
#Use nrows optional argument if you only want to read parts of the file
data=read_tsv("zcat example/GM12878_MboI_HICall_SEMA3C.tsv.gz")

#plot the whole region at 10kb resolution
plot_binned(data, resolution=10000, b1=data[,min(rbegin1)], e1=data[,max(rend2)])

#plot the whole region at 5kb resolution
plot_binned(data, resolution=5000, b1=data[,min(rbegin1)], e1=data[,max(rend2)])

#arrow plots need a category column. we can add a dummy one
data[,category:="NA"]
#plot a 20kb subset of it with base resolution (arrow plot)
plot_raw(data, b1=data[,min(rbegin1)+50000], e1=data[,min(rbegin1)+70000])

#If you want binless categories, the positions of the dangling ends and the sonication fragment size must be determined
#The following command returns a list that contains three plots and the data used to generate them
#be sure to set the proper read length, which should be
data[,median(length1)]
#if you use a Mac, use gzcat instead of zcat, or provide the path to an uncompressed file
a=examine_dataset(data, skip.fbm=T, read.len=101)
a$pdangling #select begin of dangling ends. With DpnII, one can expect 0 for dangling.L and 3 for dangling.R
a$pdiag #the distribution of sonication fragments. In this experiment, reads are mostly smaller than maxlen=900

#now we add the correct category label for each read type
data = categorize_by_new_type(data, dangling.L = c(0), dangling.R = c(3), maxlen = 900, read.len=101)

#plot the same region as before, with the new colours
plot_raw(data, b1=data[,min(rbegin1)+50000], e1=data[,min(rbegin1)+70000])

#plot a region that's further away from the diagonal
plot_raw(data, b1=data[,min(rbegin1)+120000], e1=data[,min(rbegin1)+130000],
         b2=data[,min(rbegin1)+900000], e2=data[,min(rbegin1)+910000])


