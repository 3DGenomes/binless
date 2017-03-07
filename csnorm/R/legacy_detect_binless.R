#' @include csnorm.R
NULL

#' common ops for binless detection
#' @keywords internal
estimate_binless_common = function(cs, cts, groups) {
  setkey(cts,name,id1,id2)
  cts=csnorm_predict_all(cs, cts, verbose=F)
  cts=groups[cts]
  cts=rbind(cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.close,log_mean=log_mean_cclose)],
            cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.far,log_mean=log_mean_cfar)],
            cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.up,log_mean=log_mean_cup)],
            cts[,.(groupname,ibin1,ibin2,bin1,bin2,count=contact.down,log_mean=log_mean_cdown)])
  cts[,var:=1/exp(log_mean)+1/cs@par$alpha]
  return(cts)
}

#' group counts to compute signal by cross-validated lasso
#' @keywords internal
estimate_binless_signal = function(cs, cts, groups, mat) {
  cts=estimate_binless_common(cs, cts, groups)
  cts=mat[,.(groupname,ibin1,ibin2,phi)][cts,,on=c("groupname","ibin1","ibin2")]
  cts[,phihat:=(count/exp(log_mean)-1+phi)]
  ret=cts[,.(phihat=weighted.mean(phihat,1/var),var=1/sum(1/var),ncounts=.N),
          by=c("groupname","ibin1","ibin2","bin1","bin2","phi")]
  ret[,value:=phihat/sqrt(var)]
  return(ret)
}


#' group counts to compute signal difference by cross-validated lasso
#' @keywords internal
estimate_binless_differential = function(cs, cts, groups, mat, ref) {
  cts=estimate_binless_common(cs, cts, groups)
  cts[,z:=(count/exp(log_mean)-1)]
  ret=cts[,.(zhat=weighted.mean(z,1/var),var=1/sum(1/var),ncounts=.N),
          by=c("groupname","ibin1","ibin2","bin1","bin2")]
  ret=merge(ret[groupname==ref,.(ibin1,ibin2,zhat.ref=zhat,var.ref=var)],ret[groupname!=ref],by=c("ibin1","ibin2"))
  ret=mat[groupname!=ref,.(groupname,ibin1,ibin2,delta)][ret,,on=c("groupname","ibin1","ibin2")]
  ret[,c("deltahat","var"):=list(delta + zhat-zhat.ref,var+var.ref)]
  ret[,c("value","zhat","zhat.ref","var.ref"):=list(deltahat/sqrt(var),NULL,NULL)]
  return(ret)
}
