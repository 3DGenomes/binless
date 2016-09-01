library(rstan)
library(ggplot2)
library(data.table)

bf_per_decade=5
dmin=1-0.01
dmax=150000+0.01
bpk=0.25
lambda=1
counts=cs@counts
biases=cs@biases
design=cs@design
bf_per_kb=bpk
verbose=F
iter=10000

#old
cs1=counts[,.(R=sum(contact.close+contact.down),L=sum(contact.far+contact.up)),by=pos1][,.(pos=pos1,L,R)]
cs2=counts[,.(R=sum(contact.far+contact.down),L=sum(contact.close+contact.up)),by=pos2][,.(pos=pos2,L,R)]
pos=data.table(pos=biases[,pos], key="pos")
sums=rbind(cs1,cs2)[,.(L=sum(L),R=sum(R)),keyby=pos][pos]
Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
data=list(Krow=Krow, S=biases[,.N],
          cutsites=biases[,pos], rejoined=biases[,rejoined],
          danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
          counts_sum_left=sums[,L], counts_sum_right=sums[,R],
          lambda_nu=lambda, lambda_delta=lambda)
op=optimizing(stanmodels$guess, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)

#new
cs1=counts[,.(R=sum(contact.close+contact.down),L=sum(contact.far+contact.up)),by=c("name","id1")][,.(name,id=id1,L,R)]
cs2=counts[,.(R=sum(contact.far+contact.down),L=sum(contact.close+contact.up)),by=c("name","id2")][,.(name,id=id2,L,R)]
pos=biases[,.(name,id)]
setkey(pos, name, id)
sums=rbind(cs1,cs2)[,.(L=sum(L),R=sum(R)),keyby=c("name","id")][pos]
#run optimization
Krow=round(biases[,(max(pos)-min(pos))/1000*bf_per_kb])
bbegin=c(1,biases[,.(name,row=.I)][name!=shift(name),row],biases[,.N+1])
data=list(Dsets=design[,.N], Biases=design[,uniqueN(genomic)],
          XB=array(design[,genomic],dim=design[,.N]), Krow=Krow, SD=biases[,.N], bbegin=bbegin,
          cutsitesD=biases[,pos], rejoined=biases[,rejoined],
          danglingL=biases[,dangling.L], danglingR=biases[,dangling.R],
          counts_sum_left=sums[,L], counts_sum_right=sums[,R],
          lambda_nu=lambda, lambda_delta=lambda)
op=optimizing(stanmodels$guess, data=data, as_vector=F, hessian=F, iter=iter, verbose=verbose, init=0)



#compare them
load("guess_old.RData")
sm.old=stan_model("guess_old.stan")
op.old=optimizing(sm.old, data=data.old, as_vector=F, hessian=F, iter=10000, verbose=T, init=0, init_alpha=1e-5,
                  tol_rel_grad=1e1, tol_rel_obj=1e1)

load("guess_new.RData")
sm.new=stan_model("guess_new.stan")
op.new=optimizing(sm.new, data=data.new, as_vector=F, hessian=F, iter=10000, verbose=T, init=0, init_alpha=1e-5,
                  tol_rel_grad=1e0, tol_rel_obj=1e0)

#plots
c(op.old$par$eRJ,op.new$par$eRJ[1])
c(op.old$par$eDE,op.new$par$eDE[1])
c(op.old$par$eC,op.new$par$eC[1])
#
ggplot(data.table(id=1:length(op.old$par$beta_nu), old=op.old$par$beta_nu, new=op.new$par$beta_nu[1,]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:length(op.old$par$beta_delta), old=op.old$par$beta_delta, new=op.new$par$beta_delta[1,]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:80, old=op.old$par$log_nu, new=op.new$par$log_nu[1:80]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:80, old=op.old$par$log_delta, new=op.new$par$log_delta[1:80]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:length(op.old$par$genprow), old=op.old$par$genprow, new=op.new$par$genprowD[1,]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:length(op.old$par$beta_nu_aug), old=op.old$par$beta_nu_aug, new=op.new$par$beta_nu_aug[1,]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:length(op.old$par$beta_nu_centered), old=op.old$par$beta_nu_centered, new=op.new$par$beta_nu_centered[1:length(op.old$par$beta_nu_centered)]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=1:length(op.old$par$beta_nu_diff), old=op.old$par$beta_nu_diff, new=op.new$par$beta_nu_diff[1,]))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")
#
ggplot(data.table(id=2:length(op.old$par$beta_nu), old=(op.old$par$beta_nu-shift(op.old$par$beta_nu)),
                  new=(op.new$par$beta_nu[1,])-shift(op.new$par$beta_nu[1,])))+
  geom_line(aes(id,old),colour="red")+geom_line(aes(id,new),colour="green")



