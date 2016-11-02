library(csnorm)
library(data.table)
library(ggplot2)

a=examine_dataset("/scratch/rao/mapped/GM12878_HindIII_dilution/Homo_sapiens_GM12878_HiC035_Homo_sapiens_rs_X_Peak1.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=25)
csd=read_and_prepare("/scratch/rao/mapped/GM12878_HindIII_dilution/Homo_sapiens_GM12878_HiC035_Homo_sapiens_rs_X_Peak1.tsv",
                     "data/rao_HiC035_HindIII_GM12878_Peak1_450k","GM","1",
                     enzyme="HindIII", circularize=-1, dangling.L=c(0),
                     dangling.R=c(3), maxlen=750, read.len=25, dmin=1000, save.data=T)

a=examine_dataset("/scratch/rao/mapped/GM12878_NcoI_dilution/Homo_sapiens_GM12878_HiC036_Homo_sapiens_rs_X_Peak1.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=25)
csd=read_and_prepare("/scratch/rao/mapped/GM12878_NcoI_dilution/Homo_sapiens_GM12878_HiC036_Homo_sapiens_rs_X_Peak1.tsv",
                     "data/rao_HiC036_NcoI_GM12878_Peak1_450k","GM","1",
                     enzyme="NcoI", circularize=-1, dangling.L=c(0),
                     dangling.R=c(3), maxlen=750, read.len=25, dmin=1000, save.data=T)

a=examine_dataset("/scratch/rao/mapped/GM12878_MboI_in_situ/HIC006_both_filled_map_chrX_Peak1.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
csd=read_and_prepare("/scratch/rao/mapped/GM12878_MboI_in_situ/HIC006_both_filled_map_chrX_Peak1.tsv",
                     "data/rao_HiC006_GM12878_Peak1_450k","GM","1",
                     enzyme="MboI", circularize=-1, dangling.L=c(0),
                     dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)

#here we plot the raw reads. We need to load the full csdata object, as only the one without the raw reads is returned.
#load("data/rao_HiCall_GM12878_Peak1_450k_csdata_with_data.RData")
#load("data/rao_HiC006_GM12878_Peak1_450k_csdata_with_data.RData")
#load("data/rao_HiC035_HindIII_GM12878_Peak1_450k_csdata_with_data.RData")
load("data/rao_HiC036_NcoI_GM12878_Peak1_450k_csdata_with_data.RData")
data=get_raw_reads(csd@data, csd@biases[,min(pos)], csd@biases[,max(pos)])
#plot_binned(data, resolution=10000, b1=csd@biases[,min(pos)], e1=csd@biases[,max(pos)])
plot_raw(data, b1=csd@biases[,min(pos)+10000], e1=csd@biases[,min(pos)+20000])
#ggsave(filename="images/rao_HiCall_GM12878_Peak1_450k_raw.pdf", width=10, height=8)
#ggsave(filename="images/rao_HiC006_GM12878_Peak1_450k_raw.pdf", width=10, height=8)
#ggsave(filename="images/rao_HiC035_HindIII_GM12878_Peak1_450k_raw.pdf", width=10, height=8)
ggsave(filename="images/rao_HiC036_NcoI_GM12878_Peak1_450k_raw.pdf", width=10, height=8)
