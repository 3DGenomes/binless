library(ggplot2)
library(data.table)
library(splines)
library(Matrix)
library(quadprog)

### genomic spline
#sparse spline design
splinedegree=3 #Cubic spline
Krow=100
Ncs=1000
cutsites = sort(runif(n=Ncs, min=1, max=10))
dx = 1.01*(max(cutsites)-min(cutsites))/(Krow-splinedegree)
t = min(cutsites) - dx*0.01 + dx * seq(-splinedegree,Krow-splinedegree+3)
Bsp = spline.des(cutsites, knots = t, outer.ok = T, sparse=T)$design
X = rbind(cbind(Bsp/2,Bsp/2),bdiag(Bsp,Bsp),bdiag(Bsp,Bsp))
Sm2=Diagonal(5*Ncs,x=c(rnorm(Ncs)^2*10,rnorm(2*Ncs)^2*5,rnorm(2*Ncs)^2))
ggplot(as.data.table(melt(as.matrix(X))))+geom_raster(aes(Var1,Var2,fill=value))
#
liota=cos(2*cutsites)+cutsites/10
lrho=(2+sin(cutsites))*cutsites^2/20
liota=liota-mean(liota)
lrho=lrho-mean(lrho)
y=matrix(c(10+(liota+lrho)/2,5+liota,5+lrho,1+liota,1+lrho)+rnorm(5*Ncs,sd=0.1),ncol=1)
ggplot(data.table(x=rep(cutsites,5),cat=rep(1:5,each=Ncs),y=y[,1]))+geom_point(aes(x,y))+facet_wrap(~cat)
data.table(x=rep(cutsites,5),cat=rep(1:5,each=Ncs),y=y[,1])[,mean(y),by=cat]
#
diags = list(rep(1,Krow), rep(-2,Krow))
D = bandSparse(Krow-2, Krow, k=0:2, diagonals=c(diags, diags[1]))
D = bdiag(D,D)
ggplot(as.data.table(melt(as.matrix(D))))+geom_raster(aes(Var1,Var2,fill=value))
#
A=crossprod(crossprod(Sm2,X),X)+Krow*0.1*crossprod(D)
cholA=Cholesky(A)
ggplot(as.data.table(melt(as.matrix(A))))+geom_raster(aes(Var1,Var2,fill=value))#A
ggplot(as.data.table(melt(as.matrix(A)==0)))+geom_raster(aes(Var1,Var2,fill=value))#A==0
ggplot(as.data.table(melt(as.matrix(crossprod(expand(cholA)[[1]])==0))))+geom_raster(aes(Var1,Var2,fill=value))#permuted A
ggplot(as.data.table(melt(as.matrix(expand(cholA)[[2]]==0))))+geom_raster(aes(Var1,Var2,fill=value))#factor
ggplot(as.data.table(melt(as.matrix(expand(cholA)[[1]]==0))))+geom_raster(aes(Var1,Var2,fill=value))#permutation

#
W=cbind(c(rep(0,3*Ncs),rep(1,Ncs),rep(0,Ncs)), c(rep(0,4*Ncs),rep(1,Ncs)))
ggplot(as.data.table(melt(as.matrix(W==0))))+geom_raster(aes(Var1,Var2,fill=value))
#
U=cbind(c(rep(1,Ncs),rep(0,4*Ncs)),c(rep(0,Ncs),rep(1,2*Ncs),rep(0,2*Ncs)),c(rep(0,3*Ncs),rep(1,2*Ncs)))
ggplot(as.data.table(melt(as.matrix(U==0))))+geom_raster(aes(Var1,Var2,fill=value))
#
Gamma=t(X)%*%W
Gamma=Diagonal(2*Krow) - Gamma%*%solve(t(Gamma)%*%solve(cholA,Gamma), t(solve(cholA,Gamma)))
all(abs(Gamma%*%Gamma-Gamma)<1e-5)
#
beta_y=solve(cholA, Gamma%*%t(X)%*%Sm2%*%y)
beta_U=solve(cholA, Gamma%*%t(X)%*%Sm2%*%U)
#
e=solve(t(U)%*%Sm2%*%(U-X%*%beta_U),t(U)%*%Sm2%*%(y-X%*%beta_y))
#beta = beta_y - beta_U%*%e
beta=solve(cholA, Gamma%*%t(X)%*%Sm2%*%(y-U%*%e))
t(W)%*%X%*%beta
#
mu=t(X)%*%W
mu=solve(t(mu)%*%solve(cholA,mu), t(mu)%*%solve(cholA, t(X)%*%Sm2%*%(y-U%*%e)))
#
lambda=(Krow-2)/(Krow^2*crossprod(D%*%beta))

###
result = data.table(x=rep(cutsites,5), cat=rep(1:5,each=Ncs), y=y[,1],
                    mu=as.matrix(X%*%beta)[,1],
                    mid=as.matrix(U%*%e)[,1],
                    muy=as.matrix(X%*%beta_y)[,1],
                    muU=as.matrix(X%*%beta_U%*%e)[,1])
ggplot(result)+
  geom_point(aes(x,y))+facet_wrap(~cat)+geom_line(aes(x,mu+mid),colour="green")
ggplot(result)+
  geom_point(aes(x,y))+facet_wrap(~cat)+geom_line(aes(x,mu,colour="mu"))+
  geom_line(aes(x,muU,colour="mu U"))+
  geom_line(aes(x,muy,colour="mu y"))+
  geom_line(aes(x,mid,colour="mid"))




### decay spline
splinedegree=3 #Cubic spline
Kdiag=50
dbins=10**seq(0,5,length.out=1000)
nbins=length(dbins)
dx = 1.01*(max(log(dbins))-min(log(dbins)))/(Kdiag-splinedegree)
t = min(log(dbins)) - dx*0.01 + dx * seq(-splinedegree,Kdiag-splinedegree+3)
X = spline.des(log(dbins), knots = t, outer.ok = T, sparse=T)$design
X = rbind(X,X)
U=cbind(c(rep(1,nbins),rep(0,nbins)),c(rep(0,nbins),rep(1,nbins)))
ggplot(as.data.table(melt(as.matrix(U))))+geom_raster(aes(Var1,Var2,fill=value))
Xt = cbind(U,X)
Sm2=Diagonal(2*nbins,x=c(rnorm(nbins)^2,2*rnorm(nbins)^2))
ggplot(as.data.table(melt(as.matrix(Xt))))+geom_raster(aes(Var1,Var2,fill=value))
#
w=matrix(c(1:nbins,rep(0,nbins)),ncol=1)
#
ldecay=5*exp(-(log(dbins)-5)^2/10)
y=matrix(c(5+ldecay+rnorm(nbins,sd=0.1),1+ldecay+rnorm(nbins,sd=0.1)),ncol=1)
ggplot(data.table(x=dbins,y=y[,1],dset=c(rep(c(0,1),each=nbins))))+geom_point(aes(x,y))+scale_x_log10()+facet_wrap(~dset)
data.table(x=dbins,y=y[,1],dset=c(rep(c(0,1),each=nbins)))[,mean(y),by=dset]
#
diags = list(rep(1,Kdiag), rep(-2,Kdiag))
D = bandSparse(Kdiag-2, Kdiag, k=0:2, diagonals=c(diags, diags[1]))
Dt = cbind(matrix(0,ncol=2, nrow=Kdiag-2),D)
ggplot(as.data.table(melt(as.matrix(Dt))))+geom_raster(aes(Var1,Var2,fill=value))
#
At=crossprod(crossprod(Sm2,Xt),Xt)+Kdiag*1*crossprod(Dt)
#
C=-bandSparse(Kdiag, Kdiag-1, k=c(0,-1),diagonals=list(diags[[1]],-diags[[1]]))
Ct=rbind(matrix(0,nrow=2,ncol=Kdiag), cbind(t(X)%*%w,C))
ggplot(as.data.table(melt(as.matrix(Ct)==0)))+geom_raster(aes(Var1,Var2,fill=value))
#
tmp2 = t(Xt)%*%Sm2%*%y
fit = solve.QP(At, t(Xt)%*%Sm2%*%y, -Ct, meq = 2)
betat = fit$solution
#
e=betat[1:2]
beta=betat[3:(Kdiag+2)]

###
result = data.table(x=log(dbins), y1=y[1:nbins,1], y2=y[(nbins+1):(2*nbins),1],
                    mu=as.matrix(X%*%beta)[,1],
                    mid1=rep(e[1],nbins), mid2=rep(e[2],nbins))
ggplot(result)+
  geom_point(aes(x,y1))+geom_point(aes(x,y2))+
  geom_line(aes(x,mu+mid1),colour="red")+geom_line(aes(x,mu+mid2),colour="green")
ggplot(result)+
  geom_point(aes(x,y1))+
  geom_line(aes(x,mu+mid1,colour="Xbeta+e1"))+
  geom_line(aes(x,mu,colour="Xbeta"))+
  geom_line(aes(x,mu_c,colour="Xbeta_c"))+
  geom_line(aes(x,mu_y,colour="Xbeta_y"))


