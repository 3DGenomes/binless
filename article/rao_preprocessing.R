library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)
library(scales)

a=examine_dataset("/scratch/rao/mapped/GM12878_MboI_in_situ/GM12878_MboI_HICall_Peak1.tsv",
                  skip=0L,nrows=1000000, skip.fbm=T, read.len=101)

registerDoParallel(cores=10)
chrs=c("chrX",  "chr1", "chr1", "chr7",   "chr3",  "chr4",  "chr21",      "chr21",   "chr5",    "chr12", "chr21", "chr22",  "chr1",  "chr3",     "chr3",     "chr7",      "chr21",      "chr12")
names=c("Peak1","SELP", "Talk", "SEMA3C", "FOXP1", "PARM1", "Comparison", "ADAMTS1", "ADAMTS2", "TBX3",  "Fig1C", "22qter", "Tbx19", "FOXP1ext", "FOXP1big", "SEMA3Cext", "ADAMTS1ext", "TBX3ext")
sizes=c("450k", "150k", "2M",   "1M",     "1.3M",  "600k",  "1.7M",       "2.3M",    "450k",    "1.5M",  "1M",    "1.7M",   "700k",  "2.3M",     "5M",       "1.5M",     "2.8M",       "2M")
foreach (chr=chrs, name=names, size=sizes, .errorhandling="remove") %dopar% {
  cat(chr,name,size,"\n")
  csd=read_and_prepare(paste0("zcat /scratch/rao/mapped/GM12878_MboI_in_situ/GM12878_MboI_HICall_",name,".tsv.gz"),
                       paste0("data/rao_HiCall_GM12878_",name,"_",size), "GM", "1",
                       enzyme="MboI", name=paste(name,"GM12878 all"), circularize=-1, dangling.L=c(0),
                       dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)
  csd=read_and_prepare(paste0("zcat /scratch/rao/mapped/IMR90_MboI_in_situ/IMR90_MboI_HICall_",name,".tsv.gz"),
                       paste0("data/rao_HiCall_IMR90_",name,"_",size), "IMR90", "1",
                       enzyme="MboI", name=paste(name,"IMR90 all"), circularize=-1, dangling.L=c(0),
                       dangling.R=c(3), maxlen=600, read.len=101, dmin=1000, save.data=T)
  for (run in c("odd","even")) {
    csd=read_and_prepare(paste0("zcat /scratch/rao/mapped/GM12878_MboI_in_situ/GM12878_MboI_HIC",run,"_",name,".tsv.gz"),
                         paste0("data/rao_HiC",run,"_GM12878_",name,"_",size), "GM", run,
                         enzyme="MboI", name=paste(name,"GM12878",run), circularize=-1, dangling.L=c(0),
                         dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)
  }
  for (run in c("odd","even")) {
    csd=read_and_prepare(paste0("zcat /scratch/rao/mapped/IMR90_MboI_in_situ/IMR90_MboI_HIC",run,"_",name,".tsv.gz"),
                         paste0("data/rao_HiC",run,"_IMR90_",name,"_",size), "IMR90", run,
                         enzyme="MboI", name=paste(name,"IMR90",run), circularize=-1, dangling.L=c(0),
                         dangling.R=c(3), maxlen=600, read.len=101, dmin=1000, save.data=T)
  }
}


#special case for replicate H, aligned on hg19
sub="Fig1C"
size="1M"
csd=read_and_prepare(paste0("/scratch/rao/remapping/hg19_chr21/SRR1658572_hg19_chr21_35M-36M_Fig1C.tsv"),
                     paste0("data/rao_HiC003_GM12878_hg19_",sub,"_",size), "GM", run,
                     enzyme="MboI", name=paste(sub,"replicate H"), circularize=-1, dangling.L=c(0),
                     dangling.R=c(3), maxlen=900, read.len=101, dmin=1000, save.data=T)




