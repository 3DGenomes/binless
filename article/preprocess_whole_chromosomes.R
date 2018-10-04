library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)
library(scales)

#forcato datasets
forcato=fread("~/simulations/cs_norm/article/forcato_datasets.csv", header=T)

#rao datasets
rao=fread("~/simulations/cs_norm/article/rao_datasets_annotated.csv")

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

if (F) {
  #all Aiden HindIII and NcoI samples follow this
  forcato[Dataset!="Rao",cluster_id]
  a=examine_dataset("zcat /scratch/forcato/datasets/tsv/Aiden_GM_HindIII_A_HindIII_both_map_chr1.tsv.gz | head -n10000001",
                    skip=0L,nrows=10000000, skip.fbm=T, read.len=76)
  ggplot(a$data)+geom_histogram(aes(length1),binwidth=1)
  a$data[,median(length1)] #must be 76
  a$pdangling #dangling at 1 and 4
  a$pdiag #maxlen=600
  a$pclose #dmin=1000
  
  #all Dixon_2012 and 2015 samples follow this
  forcato[Dataset!="Rao",cluster_id]
  a=examine_dataset("zcat /scratch/forcato/datasets/tsv/Dixon2012_IMR90_HindIII_B_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=36)
  ggplot(a$data)+geom_histogram(aes(length1),binwidth=1)
  a$data[,median(length1)] #must be 36
  a$pdangling #dangling at 1 and 4
  a$pdiag #maxlen=600
  a$pclose #dmin=1000

  #except Dixon2015_hESC_HindIII_B
  forcato[Dataset!="Rao",cluster_id]
  a=examine_dataset("zcat /scratch/forcato/datasets/tsv/Dixon2015_hESC_HindIII_B_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=50)
  ggplot(a$data)+geom_histogram(aes(length1),binwidth=1)
  a$data[,median(length1)] #must be 50
  a$pdangling #dangling at 1 and 4
  a$pdiag #maxlen=700
  a$pclose #dmin=1000
  
  #these Jin samples "Jin_hESC_HindIII_A"         "Jin_IMR90_HindIII_A"        "Jin_IMR90_HindIII_B" 
  forcato[Dataset!="Rao",cluster_id]
  a=examine_dataset("zcat /scratch/forcato/datasets/tsv/Jin_IMR90_HindIII_C_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=36)
  ggplot(a$data)+geom_histogram(aes(length1),binwidth=1)
  a$data[,median(length1)] #must be 36
  a$pdangling #dangling at 1 and 4
  a$pdiag #maxlen=600
  a$pclose #dmin=1000
 
  #these Jin samples "Jin_IMR90_HindIII_C"        "Jin_IMR90_HindIII_D"  "Jin_IMR90_HindIII_E"        "Jin_IMR90_HindIII_F"      
  forcato[Dataset!="Rao",cluster_id]
  a=examine_dataset("zcat /scratch/forcato/datasets/tsv/Jin_IMR90_HindIII_C_both_map_chr1.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=50)
  ggplot(a$data)+geom_histogram(aes(length1),binwidth=1)
  a$data[,median(length1)] #must be 50
  a$pdangling #dangling at 1 and 4
  a$pdiag #maxlen=700
  a$pclose #dmin=1000
  
  
  #sexton fly dataset      
  forcato[Dataset!="Rao",cluster_id]
  a=examine_dataset("zcat /scratch/forcato/datasets/tsv/Sexton_fly_DpnII_A_both_map_chr3L.tsv.gz | head -n1000001",
                    skip=0L,nrows=1000000, skip.fbm=T, read.len=36)
  ggplot(a$data)+geom_histogram(aes(length1),binwidth=1)
  a$data[,median(length1)] #must be 36
  a$pdangling #dangling at 0 and 3
  a$pdiag #maxlen=900
  a$pclose #dmin=1500
  
 

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


