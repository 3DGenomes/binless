library(data.table)

load("data/rao_HiCall_FOXP1_1.3M_bpk30_csnorm_optimized.RData")
resolution=20000

mat=get_interactions(cs, type="interactions", resolution=resolution, group="all", threshold=0.95, ref="expected")
mat=mat[,.(name,bin1,begin1,end1,bin2,begin2,end2,ncounts,observed,expected,expected.sd,
           normalized,normalized.sd,signal,signal.sd,is.significant,direction)]

#observed data and expected counts
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=log(observed)))+
  geom_raster(aes(begin2,begin1,fill=log(expected)))+facet_wrap(~name)+
  scale_fill_gradient(high="red", low="white", na.value = "white")

#normalized matrix
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=pmin(normalized,5)))+
  geom_raster(aes(begin2,begin1,fill=pmin(normalized,5)))+facet_wrap(~name)+
  scale_fill_gradient(high="red", low="white", na.value = "white")

#normalized matrix with significant interactions 
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=pmin(normalized,5)))+
  geom_raster(aes(begin2,begin1,fill=pmin(normalized,5)))+facet_wrap(~name)+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T&direction=="enriched"])+
  scale_fill_gradient(high="red", low="white", na.value = "white")

#signal matrix with significant interactions 
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=signal))+
  geom_raster(aes(begin2,begin1,fill=signal))+facet_wrap(~name)+
  geom_point(aes(begin1,begin2,colour=direction),data=mat[is.significant==T&direction=="enriched"])+
  scale_fill_gradient(high="red", low="white", na.value = "white")

#binless signal
mat=get_interactions(cs, type="binteractions", resolution=resolution, group="all", threshold=-1, ref="expected")
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=value))+
  geom_raster(aes(begin2,begin1,fill=value))+
  facet_wrap(~name)+scale_fill_gradient(high="red", low="white", na.value = "white")

#binless patches
ggplot(mat)+geom_raster(aes(begin1,begin2,fill=factor(patchno)))+
  geom_raster(aes(begin2,begin1,fill=factor(patchno)))+
  facet_wrap(~name)+guides(fill=F)+scale_colour_brewer(type="seq")
