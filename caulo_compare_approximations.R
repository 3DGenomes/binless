library(ggplot2)
library(data.table)
library(binless)
library(foreach)
library(doParallel)

setwd("/home/yannick/simulations/cs_norm")

### zoom on smaller portions
begin=500000
end=650000
load("data/caulo_NcoI_all_csdata_with_data.RData")
data=csd@data[re.closest1>=begin&re.closest1<=end&re.closest2>=begin&re.closest2<=end]
cs_data = binless:::prepare_for_sparse_cs_norm(data, both=F, circularize=-1)
csd = new("CSdata", info=csd@info,
          settings=list(circularize=-1,dmin=csd@settings$dmin,dmax=cs_data$biases[,max(pos)-min(pos)]),
          data=data, biases=cs_data$biases, counts=cs_data$counts)
save(csd, file="data/caulo_NcoI_500k_csdata_with_data.RData")
csd@data=data.table()
save(csd, file="data/caulo_NcoI_500k_csdata.RData")


#exact run
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs=run_exact(cs, bf_per_kb = 1, bf_per_decade = 10, lambdas = 10**seq(from=-3,to=0,length.out=10), ncores = 30, iter = 100000, prefix="tmp/caulo_NcoI_500k")
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_exact.RData")

#approximation run
load("data/caulo_NcoI_500k_csdata.RData")
cs=merge_cs_norm_datasets(list(csd), different.decays="none")
cs = run_gauss(cs, bf_per_kb=1, bf_per_decade=10, bins_per_bf=10, ngibbs = 20, iter=100000, init_alpha=1e-7, ncounts = 1000000, type="perf")
save(cs, file="data/caulo_NcoI_500k_csnorm_optimized_gauss_nofill_onego.RData")



#plots
dsets=c(paste0("data/caulo_NcoI_500k_csnorm_optimized_exact.RData"),
        #paste0("data/caulo_NcoI_500k_csnorm_optimized_gauss_initexact.RData"),
        #paste0("data/caulo_NcoI_500k_csnorm_optimized_exact_initgauss.RData"),
        paste0("data/caulo_NcoI_500k_csnorm_optimized_gauss.RData"))
names=c("exact",
        #"gaussie",
        #"exactig",
        "approximation")

dsets=c("data/discarded/caulo_NcoI_500k_csnorm_optimized_gauss_nofill_outer.RData",
        "data/caulo_NcoI_500k_csnorm_optimized_gauss_nofill.RData")
names=c("outer","perf")

#nu and delta
iota = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(pos=cs@biases[,pos],iota=exp(cs@par$log_iota),rho=exp(cs@par$log_rho),method=j)
}
ggplot(iota)+geom_line(aes(pos,iota,colour=method))
ggsave(filename = "images/caulo_NcoI_500k_iota_bias.pdf", width=10, height=7)
ggplot(iota)+geom_line(aes(pos,rho,colour=method))
ggsave(filename = "images/caulo_NcoI_500k_rho_bias.pdf", width=10, height=7)
#
ggplot(merge(iota[method=="exact",.(pos,iotaref=iota,rhoref=rho)],iota[method!="exact"],by="pos"))+
  geom_point(aes(iotaref,iota,colour=method))+stat_function(fun=identity)+xlab("exact")+ylab("approximation")+
  scale_x_log10(limits=c(1e-2,1e2))+scale_y_log10(limits=c(1e-2,1e2))+guides(colour=F)
ggsave(filename = "images/caulo_NcoI_500k_iota_bias_correlation.pdf", width=10, height=7)
ggplot(merge(iota[method=="exact",.(pos,iotaref=iota,rhoref=rho)],iota[method!="exact"],by="pos"))+
  geom_point(aes(rhoref,rho,colour=method))+stat_function(fun=identity)+xlab("exact")+ylab("approximation")+
  scale_x_log10(limits=c(1e-2,1e2))+scale_y_log10(limits=c(1e-2,1e2))+guides(colour=F)
ggsave(filename = "images/caulo_NcoI_500k_rho_bias_correlation.pdf", width=10, height=7)
#
cor.test(iota[method=="exact",log(iota)],iota[method=="approximation",log(iota)])
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
ggsave(filename="images/caulo_NcoI_500k_diagonal_decay.pdf", width=10, height=7)

#parameters
params = foreach(i=dsets,j=names,.combine=rbind) %do% {
  load(i)
  data.table(method=j,name=cs@design[,name],eRJ=cs@par$eRJ,eDE=cs@par$eDE,eC=cs@par$eC,alpha=cs@par$alpha,lambda_iota=cs@par$lambda_iota,
             lambda_rho=cs@par$lambda_rho,lambda_diag=cs@par$lambda_diag,value=cs@par$value)
}
params

