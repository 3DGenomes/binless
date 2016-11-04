library(ggplot2)
library(data.table)
library(csnorm)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")


### zoom on smaller portions
begin=73780000
#end=73880000 #100k
end=73840000 #60k
load("data/rao_HiCall_GM12878_Peak1_450k_csdata_with_data.RData")
data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
cs_data = csnorm:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
csd = new("CSdata", info=csd@info,
          settings=list(circularize=-1,dmin=csd@settings$dmin,dmax=cs_data$biases[,max(pos)-min(pos)]),
          data=data, biases=cs_data$biases, counts=cs_data$counts)
save(csd, file="data/rao_HiCall_GM12878_Peak1_sub60k_csdata_with_data.RData")
csd@data=data.table()
save(csd, file="data/rao_HiCall_GM12878_Peak1_sub60k_csdata.RData")


### take less reads
foreach (sub=c(1,10,50)) %do% {
  load("data/rao_HiCall_GM12878_Peak1_sub60k_csdata_with_data.RData")
  data=csd@data[sample(.N,as.integer(.N*sub/100))]
  cs_data = csnorm:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
  csd = new("CSdata", info=csd@info,
            settings=csd@settings, data=data, biases=cs_data$biases, counts=cs_data$counts)
  csd@data=data.table()
  save(csd, file=paste0("data/rao_HiCall_GM12878_Peak1_sub60k_",sub,"pc_csdata.RData"))
}


#plots
sub="sub100k"
dsets=c(paste0("data/rao_HiCall_GM12878_Peak1_",sub,"_csnorm_optimized_exact_initgauss.RData"),
        paste0("data/rao_HiCall_GM12878_Peak1_",sub,"_csnorm_optimized_gauss_initexact.RData"),
        paste0("data/rao_HiCall_GM12878_Peak1_",sub,"_csnorm_optimized_gam.RData"),
        paste0("data/rao_HiCall_GM12878_Peak1_",sub,"_csnorm_optimized_gauss.RData"))
names=c("exactig",
        "gaussie",
        "gam",
        "approximation")


dsets=c(paste0("data/rao_HiCall_GM12878_Peak1_sub60k_1pc_csnorm_optimized_exact_initgauss.RData"),
        paste0("data/rao_HiCall_GM12878_Peak1_sub60k_10pc_csnorm_optimized_exact_initgauss.RData"),
        paste0("data/rao_HiCall_GM12878_Peak1_sub60k_50pc_csnorm_optimized_exact_initgauss.RData"),
        paste0("data/rao_HiCall_GM12878_Peak1_sub60k_csnorm_optimized_exact_initgauss.RData"))
names=c("1pc",
        "10pc",
        "50pc",
        "all")


#nu and delta
iota = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],iota=exp(cs@par$log_iota),rho=exp(cs@par$log_rho),method=j)
}
ggplot(iota)+geom_line(aes(pos,iota,colour=method))
#ggsave(filename = "images/rao_HiCall_chrX_450k_iota_bias.pdf", width=10, height=7)
ggplot(iota)+geom_line(aes(pos,rho,colour=method))
#ggsave(filename = "images/rao_HiCall_chrX_450k_rho_bias.pdf", width=10, height=7)
#
ggplot(merge(iota[method=="exact",.(pos,iotaref=iota,rhoref=rho)],iota[method!="exact"],by="pos"))+
  geom_point(aes(iotaref,iota,colour=method))+stat_function(fun=identity)
ggsave(filename = "images/rao_HiCall_chrX_450k_iota_bias_correlation.pdf", width=10, height=7)
ggplot(merge(iota[method=="exact",.(pos,iotaref=iota,rhoref=rho)],iota[method!="exact"],by="pos"))+
  geom_point(aes(rhoref,rho,colour=method))+stat_function(fun=identity)
ggsave(filename = "images/rao_HiCall_chrX_450k_rho_bias_correlation.pdf", width=10, height=7)
#
cor.test(iota[method=="exact",log(iota)],iota[method=="simplified",log(iota)])
cor.test(iota[method=="exact",log(iota)],iota[method=="approximation",log(iota)])
cor.test(iota[method=="exact",log(rho)],iota[method=="simplified",log(rho)])
cor.test(iota[method=="exact",log(rho)],iota[method=="approximation",log(rho)])


#decay
decay = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  if ("decay" %in% names(cs@par$decay)) {
    cs@par$decay[,.(method=j,dist,decay)]
  } else {
    cs@par$decay[,.(method=j,dist=distance,decay=exp(kappa-cs@par$eC))]
  }
}
ggplot(decay[,.SD[sample(.N,min(.N,100000))],by=method])+
  geom_line(aes(dist,decay,colour=method))+scale_x_log10()+scale_y_log10()
#ggsave(filename="images/rao_HiCall_chrX_450k_diagonal_decay.pdf", width=10, height=7)

#parameters
params = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(method=j,eC=cs@par$eC,alpha=cs@par$alpha,lambda_iota=cs@par$lambda_iota,
             lambda_rho=cs@par$lambda_rho,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
}
params

