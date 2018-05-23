library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)
library(scales)

#forcato datasets
forcato=fread("article/forcato_datasets.csv", header=T)

#rao datasets
rao=fread("article/rao_datasets_annotated.csv")

if (F) {
  #Rao GM12878 forcato "replicate H" (HIC003 hg38)
  a=examine_dataset("zcat /scratch/forcato/rao/tsv/HIC003_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
  a$data[,median(length1)] #must be 101
  a$pdangling #dangling at 0 and 3 as expected
  a$pdiag #maxlen=900
  a$pclose #dmin=1000
  
  #Rao GM12878 second replicate (hg38)
  a=examine_dataset("zcat /scratch/forcato/rao/tsv/HIC020_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
  a$data[,median(length1)] #must be 101
  a$pdangling #dangling at 0 and 3 as expected
  a$pdiag #maxlen=900
  a$pclose #dmin=1000
  
  #Rao IMR90
  a=examine_dataset("zcat /scratch/forcato/rao/tsv/HIC050_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=101)
  a$data[,median(length1)] #must be 101
  a$pdangling #dangling at 0 and 3 as expected
  a$pdiag #maxlen=900
  a$pclose #dmin=1000
}

registerDoParallel(cores=10)
chrs=paste0("chr",c(1:22,"X"))
foreach (i=1:rao[,.N], .errorhandling="remove") %:% foreach (chr=chrs, .errorhandling="remove") %dopar% {
  cat(rao[i,sample_id],chr,"\n")
  csd=read_and_prepare(paste0("zcat /scratch/forcato/rao/tsv/",rao[i,sample_id],"_both_map_",chr,".tsv.gz"),
                       paste0("/scratch/forcato/rao/csdata/rao_",rao[i,Cell_type],"_",rao[i,RE],"_",rao[i,forcato_id],"_",rao[i,sample_id]),
                       rao[i,Cell_type], rao[i,forcato_id], enzyme=rao[i,RE],
                       name=rao[i,paste(Cell_type,RE,forcato_id)], circularize=-1,
                       dangling.L=rao[i,dangling.L], dangling.R=rao[i,dangling.R],
                       maxlen=rao[i,maxlen], read.len=rao[i,read.len], dmin=rao[i,dmin], save.data=F)
}


