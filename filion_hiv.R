library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)


setwd("/home/yannick/simulations/cs_norm")

a=examine_dataset("/scratch/fillon/E6N_hg19e6_chr2_199M-201M.tsv",skip=0,nrows=1e6,skip.fbm=F, read.len=75)
'a$data[category!="none"]
   re.closest2 re.closest2.idx re.closest1 re.closest1.idx                                      id    begin1 strand1 length1    re.up1
1:   199002357              10   199002357              10  D00733:206:CB0CFANXX:1:1215:18796:3007 199002287       0      75 199002086
2:   199003072              11   199002357              10 D00733:206:CB0CFANXX:2:2311:17438:73274 199002287       0      75 199002086
3:   199003126              12   199003072              11 D00733:209:CB0CBANXX:1:2111:15749:87306 199003073       1      54 199003072
4:   199003126              12   199003072              11 D00733:209:CB0CBANXX:2:2205:15446:26245 199003073       0      61 199003072
re.dn1    begin2 strand2 length2    re.up2    re.dn2   rbegin1     rend1   rbegin2     rend2 category
1: 199002357 199002567       0      75 199002357 199003072 199002287 199002213 199002567 199002493     test
2: 199002357 199002731       0      75 199002357 199003072 199002287 199002213 199002731 199002657     test
3: 199003126 199003188       1      75 199003126 199003655 199003073 199003126 199003188 199003262    test2
4: 199003126 199003193       0      75 199003126 199003655 199003073 199003013 199003193 199003119    test2
'
plot_raw(a$data,b1=a$data[,min(rbegin1)+2000],e1=a$data[,min(rbegin1)+4000])