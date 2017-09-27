library(csnorm)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(scales)

setwd("/home/yannick/simulations/cs_norm")

cistrans = foreach (cell=rep(c("GM12878"),18), run=sapply(1:18,function(x){sprintf("%02d",x)}), .combine=rbind) %do% {
  #read file
  data=read_tsv(paste0("/scratch/rao/mapped/",cell,"_MboI_in_situ/",cell,"_MboI_HIC0",run,"_Talk.tsv")) 
  #classify reads
  maxlen=900
  read.len=101
  dangling.L=c(0)
  dangling.R=c(3)
  #remove fbm
  data = data[!grepl("[#~]",id)]
  #remove far from diagonal and not pointing inwards
  data = data[rbegin2-rbegin1 < maxlen & strand1==1 & strand2==0 & length1 %in% read.len]
  #generate stats
  data[,is.dangling:=((rbegin1 - re.closest1) %in% dangling.L) & ((rbegin2 - re.closest2) %in% dangling.R)]
  dcast(data[,.(cell=cell,run=run,.N),by=is.dangling],cell+run~is.dangling,value.var="N")[,
                .(cell,run,ndangling=get("TRUE"),ntotal=get("TRUE")+get("FALSE"))]
}

cistrans[,dangling.pc:=100*ndangling/ntotal]
cistrans[,run:=paste0("HIC0",run)]

a=fread("/scratch/rao/mapped/GM12878_MboI_in_situ/cistrans_chr1.dat")
cistrans = merge(cistrans,a[,.(run=V2,cis=V4,trans=V6,cis.pc=V10)],by="run")
cistrans[,cell:=NULL]
ggplot(melt(cistrans,id.vars=c("run")))+geom_point(aes(run,value,colour=variable))+facet_wrap(~variable,scales="free")

ggplot(cistrans,aes(100-dangling.pc,cis.pc))+geom_point()+theme_minimal()+stat_smooth(method="lm")+
  labs(x="ligation ratio (%)", y="cis interactions (%)")
ggsave(filename="images/cistrans_vs_dangling.pdf",width=6,height=5)
cistrans[,summary(lm(cis.pc~I(100-dangling.pc)))]
cistrans[,cor.test(dangling.pc,cis.pc)]
