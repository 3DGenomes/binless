library(data.table)
library(scales)
library(ggplot2)
sfg=scale_fill_gradient(high=muted("red"), low="white", na.value = "white")
sfg2=scale_fill_gradient2(high=muted("red"), low=muted("blue"), na.value = "white")


resolution=10000 #try 10k, 20k and 50k
load(paste0("data/rao_HiCall_FOXP1_1.3M_",resolution/1000,"k_signal.RData"))



#observed data and expected counts
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+sfg+facet_wrap(~ name)

#normalized matrix
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=pmin(normalized,5)))+
  geom_raster(aes(begin2,begin1,fill=pmin(normalized,5)))+facet_wrap(~name)+sfg

#normalized matrix with significant interactions 
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=pmin(normalized,5)))+
  geom_raster(aes(begin2,begin1,fill=pmin(normalized,5)))+facet_wrap(~name)+sfg+
  geom_point(aes(begin1,begin2,colour=signal.direction),data=mat[signal.is.significant==T&signal.direction=="enriched"])

#binned signal with significant differences 
ggplot(mat[!is.na(binned.difference)])+geom_raster(aes(begin1,begin2,fill=log(binned.difference)))+
  facet_wrap(~name)+sfg2
ggplot(mat[!is.na(binned.difference)])+geom_raster(aes(begin1,begin2,fill=log(binned.difference)))+
  geom_raster(aes(begin2,begin1,fill=log(binned.difference)))+facet_wrap(~name)+sfg2+
  geom_point(aes(begin1,begin2,colour=difference.direction),data=mat[difference.is.significant==T])


#binless signal
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=binless.signal))+
  geom_raster(aes(begin2,begin1,fill=binless.signal))+
  facet_wrap(~name)+sfg

#binless differences
ggplot(mat[!is.na(binless.difference)])+facet_wrap(~name)+sfg+
  geom_raster(aes(begin1,begin2,fill=binless.difference))+
  geom_raster(aes(begin2,begin1,fill=binless.difference))
  
